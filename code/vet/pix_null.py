# Step 3b: null-calibrated joint PRF fit. For each sector fit [target PRF, G8.9 giant PRF (83"), G19.1 star (29"), const] to the in-vs-out difference image
# and to each of the ~120 random-epoch null images saved by diff.py; significance = amplitude / std of null amplitudes (accounts for the giant's own variability).
import warnings; warnings.filterwarnings('ignore')
import sys, json, numpy as np
sys.path.insert(0,'glue'); from common import *
from PRF import TESS_PRF
SRC={'target':1058588908797940864,'giant_G8.9_83as':1058589011877157632,'G19.1_29as':1058588904502228992}
req={r['sector']:r for r in __import__('pandas').read_csv('pix/out/neighbour_tests.csv').to_dict('records') if r['source']==SRC['target']}
out=[]
for s in map(int,sys.argv[1:]):
    tpf=tpf_for(s); z=np.load(f'pix/out/diff_s{s}.npz'); D,nulls,sig=z['D'],z['nulls'],z['sig']; ny,nx=D.shape
    prf=TESS_PRF(tpf.camera,tpf.ccd,s,tpf.column+nx/2,tpf.row+ny/2); g=gaia_field(tpf).set_index('Source')
    M=[prf.locate(g.loc[v].x,g.loc[v].y,(ny,nx)).ravel() for v in SRC.values()]
    A=np.vstack(M+[np.ones(D.size)]).T; W=1/sig.ravel()**2
    sol=lambda img: np.linalg.solve(A.T@(A*W[:,None]),A.T@(img.ravel()*W))
    a=sol(D); an=np.array([sol(n) for n in nulls]); sd=an.std(0,ddof=1)
    rec=dict(sector=s,n_null=len(nulls),A_req_target=req[s]['A_req'] if s in req else None)
    for i,k in enumerate(SRC): rec[k]=dict(amp=float(a[i]),null_sd=float(sd[i]),z=float(a[i]/sd[i]),frac_null_ge=float(np.mean(an[:,i]>=a[i])))
    out.append(rec); print(s,{k:(round(rec[k]['amp'],1),round(rec[k]['null_sd'],1),round(rec[k]['z'],1)) for k in SRC},'A_req_target',round(rec['A_req_target'],2),flush=True)
json.dump(out,open('evidence/03b_pixel_null_joint.json','w'),indent=1)
