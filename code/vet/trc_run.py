# TRICERATOPS FPP/NFPP, no neighbours removed. Adapted from v4/vet_264151965/glue/trc_run(_molusc).py. usage: trc_run.py RUN_INDEX [molusc]   (RUN_INDEX=-1 builds the target only)
# Light curve: native-cadence PDCSAP windows (evidence/01a_windows.npy, SPOC 2-min S14/S21 + TESS-SPOC FFI), folded at the free-fit P/T0, 10-min bins, +-0.4 d; exptime 10 min.
# Apertures: the pipeline aperture (bit 2) from each sector's light-curve file.
import sys, json, pickle, os, time, glob, numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
sys.path.insert(0,'/Users/advaitpaliwal/Scratch/planet-hunt/hunt-ffi/trc'); sys.path.insert(0,'/home/advaitpaliwal/planet-hunt/trc_264151965/glue'); import gaia_patch
import triceratops.triceratops as tr
tic=149928372; k=int(sys.argv[1]); MOL=len(sys.argv)>2 and sys.argv[2]=='molusc'; od='trc_molusc/out' if MOL else 'trc/out'
c=json.load(open('evidence/01b_fit_free.json')); P,T0=c['P'][1],c['T0'][1]; T14=c['T14_h'][1]/24
A=np.load('evidence/01a_windows.npy'); t,f=A[:,0],A[:,1]; secs=[14,21,41,47,48,74,75]
pk='trc/out/target.pkl'
if not os.path.exists(pk):
    tg=tr.target(ID=tic,sectors=np.array(secs)); aps=[]
    for i,s in enumerate(secs):
        fn=glob.glob(f'lcfits/*_s{s}.fits')
        with fits.open(fn[0]) as h:
            a=h['APERTURE'].data.astype(int); c0=h['APERTURE'].header['CRVAL1P']; r0=h['APERTURE'].header['CRVAL2P']
        rows,cols=np.nonzero(np.bitwise_and(a,2)); aps.append(np.column_stack([cols+c0,rows+r0]))
    tg._aps=aps; pickle.dump(tg,open(pk,'wb'))
    tg.stars.to_csv('trc/out/stars.csv',index=False); print('target built',len(tg.stars),'stars',flush=True)
tg=pickle.load(open(pk,'rb'))
ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.4
edges=np.arange(-0.4,0.4+1e-9,10/1440); idx=np.digitize(ph[w],edges)
tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
oot=np.abs(tb)>T14*0.75; err=float(np.std(fb[oot])); depth=float(np.median(fb[oot])-fb[np.abs(tb)<T14/2*0.8].mean())
tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps)
if k<0:
    tg.stars.to_csv('trc/out/stars_with_tdepth.csv',index=False); print('depth',depth,'err',err,'nbins',len(tb),'tdepth>0',int((tg.stars.tdepth>0).sum())); print(tg.stars.head(15).to_string()); sys.exit()
MF=os.path.abspath('molusc_out/molusc_for_trc_N1e6.csv') if MOL else None
np.random.seed(1000+k); t0=time.time()
kw=dict(molusc_file=MF) if MOL else {}
tg.calc_probs(time=tb,flux_0=fb,flux_err_0=err,P_orb=P,contrast_curve_file=None,filt='TESS',N=1000000,parallel=True,exptime=10/1440,nsamples=20,verbose=0,**kw)
tg.probs.to_csv(f'{od}/probs_{k}.csv',index=False)
json.dump(dict(tic=tic,run=k,molusc=MOL,P=P,T0=T0,depth_ppm=depth*1e6,err_ppm=err*1e6,nbins=len(tb),n_stars=len(tg.stars),n_tdepth_gt0=int((tg.stars.tdepth>0).sum()),FPP=float(tg.FPP),NFPP=float(tg.NFPP),sec=time.time()-t0),open(f'{od}/res_{k}.json','w'))
print(k,tg.FPP,tg.NFPP,round(time.time()-t0),flush=True)
