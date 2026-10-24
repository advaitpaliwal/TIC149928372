# Step 1c: plots + per-epoch depths + odd/even + phase-0.5 on native-cadence windows, using the free-fit ephemeris; compares free vs density-prior fits.
import numpy as np, json, batman, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A=np.load('evidence/01a_windows.npy'); t,f,e,ex,sec,ep=A.T
F=json.load(open('evidence/01b_fit_free.json')); R=json.load(open('evidence/01b_fit_prior.json'))
P,T0=F['P'][1],F['T0'][1]; T14=F['T14_h'][1]/24
ph=((t-T0)/P+0.5)%1-0.5; dt=ph*P; n=np.round((t-T0)/P).astype(int)
def box(m):
    a=m&(np.abs(dt)<0.4*T14); o=m&(np.abs(dt)>0.75*T14)
    w=1/e**2; d=(np.sum(w[o]*f[o])/w[o].sum()-np.sum(w[a]*f[a])/w[a].sum()); s=np.sqrt(1/w[a].sum()+1/w[o].sum()); return float(d*1e6),float(s*1e6),int(a.sum())
per=[]
for k in np.unique(n):
    d,s,na=box(n==k); per.append(dict(epoch=int(k),sector=int(sec[n==k][0]),depth_ppm=d,err_ppm=s,n_in=na))
D=np.array([p['depth_ppm'] for p in per]); S=np.array([p['err_ppm'] for p in per]); wm=np.sum(D/S**2)/np.sum(1/S**2)
odd=box(n%2==1); even=box(n%2==0)
out=dict(P=P,T0=T0,T14_h=T14*24,per_epoch=per,weighted_mean_ppm=float(wm),chi2=float(np.sum(((D-wm)/S)**2)),dof=len(D)-1,
  odd=odd,even=even,odd_even_sigma=abs(odd[0]-even[0])/np.hypot(odd[1],even[1]),
  sap_note='windows were detrended with a quadratic outside +-0.15 d, so phase 0.5 is checked in s1_secondary.py on the search npz')
json.dump(out,open('evidence/01c_epochs_oddeven.json','w'),indent=1); print(json.dumps({k:v for k,v in out.items() if k!='per_epoch'},indent=1))
def mod(fit,x):
    b=fit['best']; pr=batman.TransitParams(); pr.t0=0;pr.per=P;pr.rp=b['rp_rs'];pr.a=b['a_rs'];pr.inc=np.degrees(np.arccos(b['b']/b['a_rs']));pr.ecc=0;pr.w=90;pr.u=fit['LD_fixed'];pr.limb_dark='quadratic'
    return batman.TransitModel(pr,x).light_curve(pr)
xs=np.linspace(-0.3,0.3,1500)
fig,ax=plt.subplots(1,2,figsize=(13,4.6))
for a_,sel,ti in [(ax[0],np.ones_like(t,bool),'all 7 sectors'),(ax[1],(sec==14)|(sec==21),'SPOC 2-min only (S14, S21)')]:
    be=np.arange(-0.3,0.3001,10/1440); idx=np.digitize(dt[sel],be); x=dt[sel]; y=f[sel]
    a_.plot(x*24,(y-1)*1e6,',',color='0.85')
    bc=[(be[j-1]+be[j])/2 for j in range(1,len(be)) if (idx==j).sum()>3]; bm=[np.mean(y[idx==j]) for j in range(1,len(be)) if (idx==j).sum()>3]; bs=[np.std(y[idx==j])/np.sqrt((idx==j).sum()) for j in range(1,len(be)) if (idx==j).sum()>3]
    a_.errorbar(np.array(bc)*24,(np.array(bm)-1)*1e6,np.array(bs)*1e6,fmt='ko',ms=3)
    a_.plot(xs*24,(mod(F,xs)-1)*1e6,'C3',label=f"free: b {F['b'][1]:.2f}, rho {F['rho_star_gcc'][1]:.2f} g/cc")
    a_.plot(xs*24,(mod(R,xs)-1)*1e6,'C0--',label=f"TIC density prior: b {R['b'][1]:.2f}, rho {R['rho_star_gcc'][1]:.2f}")
    a_.set_ylim(-2500,1500); a_.set_xlim(-7,7); a_.legend(fontsize=8); a_.set_xlabel('hours from mid-transit'); a_.set_ylabel('ppm'); a_.set_title(f'TIC 149928372, P {P:.5f} d, {ti} (10-min bins)')
plt.tight_layout(); plt.savefig('plots/01_transit_fit.png',dpi=110)
fig,a_=plt.subplots(figsize=(10,4)); a_.errorbar([p['epoch'] for p in per],D,S,fmt='o'); a_.axhline(wm,ls=':')
for p in per: a_.annotate(f"S{p['sector']}",(p['epoch'],p['depth_ppm']+p['err_ppm']+40),fontsize=7,ha='center')
a_.set_xlabel('epoch'); a_.set_ylabel('depth ppm (box, central 80% of T14)'); a_.set_title(f"per-transit depth, chi2 {out['chi2']:.1f}/{out['dof']}"); plt.tight_layout(); plt.savefig('plots/01_per_epoch.png',dpi=110)
fig,axs=plt.subplots(3,5,figsize=(16,8)); axs=axs.ravel()
for i,k in enumerate(np.unique(n)[:15]):
    m=n==k; axs[i].plot(dt[m]*24,(f[m]-1)*1e6,'.',ms=1.5,color='0.6'); be=np.arange(-0.6,0.6001,20/1440); idx=np.digitize(dt[m],be)
    axs[i].plot([(be[j-1]+be[j])/2*24 for j in range(1,len(be)) if (idx==j).sum()>1],[(np.mean(f[m][idx==j])-1)*1e6 for j in range(1,len(be)) if (idx==j).sum()>1],'k.-')
    axs[i].set_ylim(-4000,3000); axs[i].axvspan(-T14*12,T14*12,color='C1',alpha=.15); axs[i].set_title(f'E{k} S{int(sec[m][0])}',fontsize=8)
for j in range(i+1,15): axs[j].axis('off')
plt.tight_layout(); plt.savefig('plots/01_every_epoch.png',dpi=90)
