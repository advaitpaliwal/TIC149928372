# Step 4: per-event depth in raw SAP vs PDCSAP (SPOC 2-min S14,S21; TESS-SPOC FFI S41-75); box = central 80% of fitted T14; local linear baseline per event, no detrending.
import numpy as np, glob, json, re, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
_c=json.load(open('evidence/01b_fit_free.json')); P,T0=_c['P'][1],_c['T0'][1]; W=0.8*_c['T14_h'][1]/24; b=dict(depth_ppm=(1-0)*0)
def events(t,f):
    ok=np.isfinite(t)&np.isfinite(f)&(f>0); t,f=t[ok],f[ok]/np.nanmedian(f[ok]); ev=[]
    for n in range(int(np.floor((t.min()-T0)/P)),int(np.ceil((t.max()-T0)/P))+1):
        dt=t-(T0+n*P); a=np.abs(dt)<W/2; bb=(np.abs(dt)>0.75*W)&(np.abs(dt)<0.75*W+0.6)
        if a.sum()<1 or bb.sum()<10 or (dt[bb]<0).sum()<3 or (dt[bb]>0).sum()<3: continue
        c=np.polyfit(dt[bb],f[bb],1); r=f[bb]-np.polyval(c,dt[bb]); s=1.4826*np.median(np.abs(r))
        ev.append(dict(epoch=n,depth_ppm=float((np.polyval(c,dt[a])-f[a]).mean()*1e6),err_ppm=float(s/np.sqrt(a.sum())*1e6),n_in=int(a.sum())))
    return ev
def comb(ev):
    if not ev: return None
    d=np.array([e['depth_ppm'] for e in ev]); e=np.array([e['err_ppm'] for e in ev]); w=1/e**2
    return dict(depth_ppm=float((d*w).sum()/w.sum()),err_ppm=float(1/np.sqrt(w.sum())),n_ev=len(ev))
out={}; allr={'raw':[], 'raw_crowd':[], 'corr':[]}
for fn in sorted(glob.glob('lcfits/*.fits'),key=lambda x:int(re.search(r'_s(\d+)',x).group(1))):
    h=fits.open(fn); d=h[1].data; s=int(re.search(r'_s(\d+)',fn).group(1)); spoc='PDCSAP_FLUX' in d.columns.names
    q=(d['QUALITY']==0); t=d['TIME'][q]
    raw=events(t,np.array(d['SAP_FLUX'][q],float)); cor=events(t,np.array(d['PDCSAP_FLUX' if spoc else 'DET_FLUX'][q],float))
    rec=dict(product='TESS-SPOC' if spoc else 'QLP',raw=comb(raw),corrected=comb(cor),raw_events=raw,corrected_events=cor)
    if spoc and rec['raw']:
        c=h[1].header['CROWDSAP']; rec['CROWDSAP']=c; rec['raw_div_crowdsap']=dict(depth_ppm=rec['raw']['depth_ppm']/c,err_ppm=rec['raw']['err_ppm']/c)
    out[f'S{s}']=rec; h.close()
    print(s,rec['product'],{k:(round(v['depth_ppm']),round(v['err_ppm'])) for k,v in rec.items() if isinstance(v,dict) and 'depth_ppm' in v})
def tot(key,prod=None):
    vals=[(v[key]['depth_ppm'],v[key]['err_ppm']) for v in out.values() if v.get(key) and (prod is None or v['product']==prod)]
    d=np.array(vals); w=1/d[:,1]**2; return dict(depth_ppm=float((d[:,0]*w).sum()/w.sum()),err_ppm=float(1/np.sqrt(w.sum())),n_sectors=len(d))
out['combined']=dict(spoc_raw=tot('raw','TESS-SPOC'),spoc_raw_div_crowdsap=tot('raw_div_crowdsap'),spoc_pdcsap=tot('corrected','TESS-SPOC'))
ks0=[k for k in out if k not in('combined','S75')]
out['combined_without_S75']={key:(lambda v:dict(depth_ppm=float((v[:,0]/v[:,1]**2).sum()/(1/v[:,1]**2).sum()),err_ppm=float(1/np.sqrt((1/v[:,1]**2).sum()))))(np.array([(out[k][key]['depth_ppm'],out[k][key]['err_ppm']) for k in ks0])) for key in ['raw','raw_div_crowdsap','corrected']}  # S75: its one event sits at a data gap (plots/01_every_epoch.png)
print(json.dumps(out['combined'],indent=1))
json.dump(out,open('evidence/04_sap_vs_pdcsap.json','w'),indent=1)
ks=[k for k in out if not k.startswith('combined')]; x=np.arange(len(ks))
fig,a=plt.subplots(figsize=(11,4))
for off,key,lab in [(-0.2,'raw','raw SAP'),(0,'raw_div_crowdsap','SAP / CROWDSAP'),(0.2,'corrected','PDCSAP or QLP DET')]:
    v=[out[k].get(key) for k in ks]; a.errorbar(x[[i for i,z in enumerate(v) if z]]+off,[z['depth_ppm'] for z in v if z],[z['err_ppm'] for z in v if z],fmt='o',label=lab)
a.set_xticks(x); a.set_xticklabels(ks); a.axhline(0,color='0.6'); a.legend(); a.set_ylabel('depth ppm'); a.set_title('per-sector transit depth, raw vs corrected flux')
plt.tight_layout(); plt.savefig('plots/04_sap_vs_pdcsap.png',dpi=110)
