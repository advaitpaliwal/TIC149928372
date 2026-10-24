# Step 1a: native-cadence PDCSAP windows (+-0.6 d) around each predicted transit, local quadratic baseline fit outside +-0.15 d, saved for the transit fit.
import numpy as np, glob, json, re, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
r=json.load(open('evidence/149928372.json'))['peaks'][0]; P,T0=r['P'],r['T0']
rows=[]
for fn in sorted(glob.glob('lcfits/*.fits')):
    s=int(re.search(r'_s(\d+)',fn).group(1))
    with fits.open(fn) as h:
        d=h[1].data; q=(d['QUALITY']==0)&np.isfinite(d['PDCSAP_FLUX'])&np.isfinite(d['TIME'])
        t=np.array(d['TIME'][q],float); f=np.array(d['PDCSAP_FLUX'][q],float); exp=float(np.median(np.diff(t)))
    for n in range(int(np.floor((t.min()-T0)/P)),int(np.ceil((t.max()-T0)/P))+1):
        dt=t-(T0+n*P); w=np.abs(dt)<0.6; o=w&(np.abs(dt)>0.15)
        if (np.abs(dt[w])<0.06).sum()<3 or (dt[o]<0).sum()<10 or (dt[o]>0).sum()<10: continue
        c=np.polyfit(dt[o],f[o],2); ff=f[w]/np.polyval(c,dt[w])
        sig=1.4826*np.median(np.abs(ff[np.abs(dt[w])>0.15]-1))
        for a,b in zip(t[w],ff): rows.append((a,b,sig,exp,s,n))
A=np.array(rows); np.save('evidence/01a_windows.npy',A)
out={}
for s in np.unique(A[:,4]):
    m=A[:,4]==s; out[int(s)]=dict(cadence_min=float(A[m,3][0]*1440),epochs=sorted(set(int(x) for x in A[m,5])),noise_ppm=float(np.median(A[m,2])*1e6))
json.dump(dict(P_search=P,T0_search=T0,sectors=out),open('evidence/01a_windows.json','w'),indent=1); print(json.dumps(out,indent=1))
