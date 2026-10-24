# Step 1d: (a) grazing test: best chi2 with b fixed at 0.85/0.9/0.95/1.0/1.05 (rp, a/Rs, T0 free; P fixed) vs free best;
#          (b) secondary eclipse: depth at phase 0.5 and over a phase grid, on native PDCSAP with the same per-window quadratic baseline as s1_prep.py.
import numpy as np, json, glob, re, batman, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from astropy.io import fits
A=np.load('evidence/01a_windows.npy'); t,f,e,ex=A[:,0],A[:,1],A[:,2],A[:,3]
F=json.load(open('evidence/01b_fit_free.json')); P=F['P'][1]; T0=F['T0'][1]; U=F['LD_fixed']
groups=[]
for c in np.unique(np.round(ex*1440,1)):
    m=np.round(ex*1440,1)==c; pr=batman.TransitParams(); pr.t0=T0;pr.per=P;pr.rp=.036;pr.a=36;pr.inc=89.9;pr.ecc=0;pr.w=90;pr.u=U;pr.limb_dark='quadratic'
    groups.append((m,pr,batman.TransitModel(pr,t[m],supersample_factor=1 if c<2.5 else int(np.ceil(c)),exp_time=c/1440)))
def chi2(tc,rp,b,a):
    if not(0.005<rp<0.6 and 1.5<a<300 and b<a): return 1e12
    mo=np.ones_like(t)
    for m,pr,TM in groups: pr.t0=tc;pr.rp=rp;pr.a=a;pr.inc=np.degrees(np.arccos(b/a)); mo[m]=TM.light_curve(pr)
    return float(np.sum(((f-mo)/e)**2))
bb=F['best']; c_free=min(minimize(lambda p:chi2(p[0],p[1],p[2],p[3]),[bb['T0'],bb['rp_rs'],b0,a0],method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-7,fatol=1e-4)).fun for b0,a0 in [(0.05,38),(0.4,35)])
G=6.674e-8; res=[]
for bf in [0.85,0.9,0.95,1.0,1.05]:
    best=None
    for rp0 in [0.04,0.06,0.1]:
        for a0 in [10,16,25]:
            r=minimize(lambda p:chi2(p[0],p[1],bf,p[2]),[T0,rp0,a0],method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-7,fatol=1e-4))
            if best is None or r.fun<best.fun: best=r
    a=best.x[2]; res.append(dict(b=bf,rp_rs=float(best.x[1]),a_rs=float(a),Rp_Re=float(best.x[1]*0.66*109.1),rho_star_gcc=float(3*np.pi*a**3/(G*(P*86400)**2)),chi2=float(best.fun),dchi2_vs_free=float(best.fun-c_free)))
    print(res[-1],flush=True)
# secondary: windows around each phase in a grid, box of T14
T14=F['T14_h'][1]/24; rows=[]
for fn in sorted(glob.glob('lcfits/*.fits')):
    with fits.open(fn) as h:
        d=h[1].data; q=(d['QUALITY']==0)&np.isfinite(d['PDCSAP_FLUX'])
        tt=np.array(d['TIME'][q],float); ff=np.array(d['PDCSAP_FLUX'][q],float); ff/=np.nanmedian(ff)
    rows.append((tt,ff))
def depth_at(phase):
    D=[];W=[]
    for tt,ff in rows:
        for n in range(int(np.floor((tt.min()-T0)/P))-1,int(np.ceil((tt.max()-T0)/P))+1):
            tc=T0+(n+phase)*P; dt=tt-tc; w=np.abs(dt)<0.6; o=w&(np.abs(dt)>0.75*T14); i=np.abs(dt)<0.4*T14
            if i.sum()<3 or (dt[o]<0).sum()<10 or (dt[o]>0).sum()<10: continue
            c=np.polyfit(dt[o],ff[o],2); r=ff[o]/np.polyval(c,dt[o])-1; s=1.4826*np.median(np.abs(r))
            D.append(1-np.mean(ff[i]/np.polyval(c,dt[i]))); W.append(i.sum()/s**2)
    D=np.array(D);W=np.array(W); return float(np.sum(D*W)/W.sum()*1e6), float(1/np.sqrt(W.sum())*1e6), len(D)
ph05=depth_at(0.5); grid=[(float(p),)+depth_at(p) for p in np.arange(0.05,0.951,0.01)]
gs=np.array([g[1]/g[2] for g in grid]); ib=int(np.argmax(gs))
out=dict(free_chi2=c_free,grazing_fixed_b=res,secondary_phase05=dict(depth_ppm=ph05[0],err_ppm=ph05[1],n_events=ph05[2],sigma=ph05[0]/ph05[1]),
  phase_grid_max=dict(phase=grid[ib][0],depth_ppm=grid[ib][1],err_ppm=grid[ib][2],sigma=float(gs[ib])),phase_grid_scatter_sigma=float(np.std(gs)),grid=grid)
json.dump(out,open('evidence/01d_shape_secondary.json','w'),indent=1); print(json.dumps({k:v for k,v in out.items() if k!='grid'},indent=1))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,a=plt.subplots(figsize=(10,3.8)); g=np.array(grid); a.errorbar(g[:,0],g[:,1],g[:,2],fmt='o',ms=3); a.axhline(0,color='0.6'); a.axvline(0.5,ls=':')
a.set_xlabel('orbital phase'); a.set_ylabel('box depth (ppm), T14 box'); a.set_title(f'search for a secondary eclipse: phase 0.5 = {ph05[0]:.0f} +- {ph05[1]:.0f} ppm (primary ~1470 ppm)'); plt.tight_layout(); plt.savefig('plots/01_secondary_phases.png',dpi=110)
