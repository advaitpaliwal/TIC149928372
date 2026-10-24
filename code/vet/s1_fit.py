# Step 1b: batman + emcee fit on native-cadence PDCSAP windows. usage: s1_fit.py MODE  (free | prior | free2m)
#  free   : rp, b, a/Rs (log-uniform), T0, P free; all 7 sectors.
#  prior  : same, Gaussian prior on stellar density from TIC 8.2 (rho 2.1562 +- 0.8084 solar = 3.04 +- 1.14 g/cc).
#  free2m : free fit, SPOC 2-min sectors only (S14, S21), P fixed to the all-data free fit if present else search P.
# Quadratic LD fixed u1 0.44 u2 0.20 (TESS band, Teff 4000 K, logg 4.6; inferred approximate Claret 2017 value, not looked up).
import numpy as np, sys, json, batman, emcee, os
mode=sys.argv[1]; A=np.load('evidence/01a_windows.npy'); t,f,e,ex,sec,ep=A.T
if mode=='free2m': k=(sec==14)|(sec==21); t,f,e,ex,sec,ep=t[k],f[k],e[k],ex[k],sec[k],ep[k]
s0=json.load(open('evidence/01a_windows.json')); P0,T00=s0['P_search'],s0['T0_search']
NREF=51; Tref0=T00+NREF*P0
U=[0.44,0.20]; RHO_SUN=1.41; rhoT,rhoTe=2.1562*RHO_SUN,0.8084*RHO_SUN; G=6.674e-8
groups=[]
for c in np.unique(np.round(ex*1440,1)):
    m=np.round(ex*1440,1)==c; pr=batman.TransitParams(); pr.t0=0;pr.per=P0;pr.rp=.03;pr.a=30;pr.inc=89.9;pr.ecc=0;pr.w=90;pr.u=U;pr.limb_dark='quadratic'
    ss=1 if c<2.5 else int(np.ceil(c/1.0)); groups.append((m,pr,batman.TransitModel(pr,t[m],supersample_factor=ss,exp_time=c/1440)))
fixP=None
if mode=='free2m' and os.path.exists('evidence/01b_fit_free.json'): fixP=json.load(open('evidence/01b_fit_free.json'))['P'][1]
def model(p):
    tc,P,rp,b,la=p; a=np.exp(la); mo=np.ones_like(t)
    for m,pr,TM in groups:
        pr.t0=tc;pr.per=P;pr.rp=rp;pr.a=a;pr.inc=np.degrees(np.arccos(b/a)); mo[m]=TM.light_curve(pr)
    return mo
def lnp(p):
    tc,P,rp,b,la=p
    if fixP is not None: P=fixP
    if not(abs(tc-Tref0)<0.1 and abs(P-P0)<0.01 and 0.005<rp<0.6 and 0<=b<1+rp and np.log(2)<la<np.log(300)): return -np.inf
    lp=0.
    if mode=='prior':
        rho=3*np.pi*np.exp(la)**3/(G*(P*86400)**2); lp=-0.5*((rho-rhoT)/rhoTe)**2
    return lp-0.5*np.sum(((f-model([tc,P,rp,b,la]))/e)**2)
nw=40; rng=np.random.default_rng(3)
p0=np.array([Tref0,P0,0.03,0.4,np.log(30)])+np.array([2e-3,2e-5,2e-3,0.25,0.15])*rng.standard_normal((nw,5)); p0[:,3]=np.abs(p0[:,3])
S=emcee.EnsembleSampler(nw,5,lnp); S.run_mcmc(p0,8000,progress=False); ch=S.get_chain(discard=3000,thin=10,flat=True)
if fixP is not None: ch[:,1]=fixP
lp=S.get_log_prob(discard=3000,thin=10,flat=True); best=ch[np.argmax(lp)]
a=np.exp(ch[:,4]); rho=3*np.pi*a**3/(G*(ch[:,1]*86400)**2); inc=np.arccos(ch[:,3]/a)
T14=ch[:,1]/np.pi*np.arcsin(np.clip(np.sqrt(np.clip((1+ch[:,2])**2-ch[:,3]**2,0,None))/a/np.sin(inc),0,1))*24
T23=ch[:,1]/np.pi*np.arcsin(np.clip(np.sqrt(np.clip((1-ch[:,2])**2-ch[:,3]**2,0,None))/a/np.sin(inc),0,1))*24
q=lambda v:[float(np.percentile(v,16)),float(np.percentile(v,50)),float(np.percentile(v,84))]
aT=(rhoT*G*(P0*86400)**2/(3*np.pi))**(1/3)
chi2best=float(np.sum(((f-model(best))/e)**2))
out=dict(mode=mode,n_pts=int(len(t)),sectors=sorted(set(int(x) for x in sec)),LD_fixed=U,Tref_epoch=NREF,T0=q(ch[:,0]),P=q(ch[:,1]),rp_rs=q(ch[:,2]),b=q(ch[:,3]),a_rs=q(a),
  rho_star_gcc=q(rho),rho_TIC_gcc=[rhoT,rhoTe],a_rs_TIC=aT,T14_h=q(T14),T23_over_T14=q(np.where(T14>0,T23/T14,0)),
  Rp_Re=q(ch[:,2]*0.66*109.1),Rp_Re_incl_Rstar_err=q(ch[:,2]*rng.normal(0.66,0.067,len(ch))*109.1),
  frac_grazing_b_gt_1_minus_rp=float(np.mean(ch[:,3]>1-ch[:,2])),frac_b_gt_0p8=float(np.mean(ch[:,3]>0.8)),frac_rp_gt_3Re=float(np.mean(ch[:,2]*0.66*109.1>3)),
  best=dict(zip(['T0','P','rp_rs','b','a_rs'],[float(best[0]),float(best[1]),float(best[2]),float(best[3]),float(np.exp(best[4]))])),chi2_best=chi2best,
  acceptance=float(np.mean(S.acceptance_fraction)))
try: out['autocorr']=S.get_autocorr_time(discard=3000,quiet=True).tolist()
except Exception as z: out['autocorr']=str(z)
np.save(f'evidence/01b_chain_{mode}.npy',ch)
json.dump(out,open(f'evidence/01b_fit_{mode}.json','w'),indent=1); print(json.dumps(out,indent=1))
