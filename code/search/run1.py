# usage: LC_ONLY=1 run1.py TIC [pool.csv|control.csv]   (run from the lane folder; LC_ONLY is the only mode)
# hz-near's mixed-product reader: one product per sector from pool.csv 'sectors' ([sector, provenance, exptime, MAST uri];
# SPOC 120-s PDCSAP, else TESS-SPOC PDCSAP, else QLP KSPSAP), quality 'hard', wotan biweight 0.75 d, 3-sigma upper clip,
# 2-min data pre-binned to 10 min (the search bins to 30 min). Saves lc/<tic>.npz (t,f,s). Each FITS is deleted right after reading.
# Exits non-zero with no npz if fewer than 3 sectors are usable.
import sys,os,json,warnings,subprocess,shutil,numpy as np,pandas as pd
warnings.filterwarnings('ignore')
import lightkurve as lk
from wotan import flatten
tic=int(sys.argv[1]);src=sys.argv[2] if len(sys.argv)>2 else 'pool.csv'
if not os.path.exists(src) or tic not in set(pd.read_csv(src,usecols=['ID']).ID):src='control.csv'
row=pd.read_csv(src).set_index('ID').loc[tic]
d=os.path.join('lc',str(tic));os.makedirs(d,exist_ok=True)
def get(uri,fn):
    r=subprocess.run(['curl','-sfL','-m','120','--retry','4','--retry-max-time','300','-o',fn,'https://mast.stsci.edu/api/v0.1/Download/file?uri='+uri])
    if r.returncode!=0:
        if os.path.exists(fn):os.remove(fn)
        return None
    return fn
T=[];F=[];SID=[]
try:
    for s,prov,exp,uri in json.loads(row.sectors):
        fn=get(uri,f'{d}/s{s}.fits')
        if not fn:continue
        lc=None
        for col in {'SPOC':['pdcsap_flux'],'TESS-SPOC':['pdcsap_flux'],'QLP':['kspsap_flux','det_flux','sap_flux']}[prov]:
            try:lc=lk.read(fn,flux_column=col,quality_bitmask='hard');break
            except Exception:lc=None
        os.remove(fn)
        if lc is None:continue
        lc=lc.remove_nans().normalize()
        if len(lc)<100:continue
        tt=lc.time.value;ff=np.asarray(lc.flux.value,float)
        fl=flatten(tt,ff,method='biweight',window_length=0.75,break_tolerance=0.5,edge_cutoff=0.0)
        m=np.isfinite(fl);tt,fl=tt[m],fl[m]
        med=np.median(fl);mad=1.4826*np.median(np.abs(fl-med));m=fl<med+3*mad;tt,fl=tt[m],fl[m]
        if np.median(np.diff(tt))<9/1440:  # 2-min: bin to 10 min
            b=np.floor((tt-tt[0])/(10/1440)).astype(int);n=np.bincount(b);k=n>0
            tt=np.bincount(b,tt)[k]/n[k];fl=np.bincount(b,fl)[k]/n[k]
        T.append(tt);F.append(fl);SID.append(np.full(len(tt),s))
finally:
    shutil.rmtree(d,ignore_errors=True)
if len(T)<3:sys.exit(f'only {len(T)} usable sectors')
t=np.concatenate(T);f=np.concatenate(F);sid=np.concatenate(SID);o=np.argsort(t)
np.savez_compressed(f'lc/{tic}.npz',t=t[o],f=f[o],s=sid[o])
