import json,glob,pandas as pd,numpy as np
out={}
for d in ['trc/out','trc_molusc/out']:
    R=[json.load(open(f)) for f in sorted(glob.glob(d+'/res_*.json'))]
    if not R: continue
    F=np.array([r['FPP'] for r in R]); N=np.array([r['NFPP'] for r in R])
    P=pd.concat([pd.read_csv(f) for f in glob.glob(d+'/probs_*.csv')])
    sc=P.groupby('scenario').prob.sum()/len(R); top=P.groupby(['ID','scenario']).prob.sum().div(len(R)).sort_values(ascending=False).head(8)
    out[d]=dict(n_runs=len(R),FPP=F.round(4).tolist(),NFPP=N.round(4).tolist(),FPP_mean=float(F.mean()),FPP_sd=float(F.std(ddof=1)),NFPP_mean=float(N.mean()),NFPP_sd=float(N.std(ddof=1)),molusc=R[0]['molusc'],depth_ppm=R[0]['depth_ppm'],n_stars=R[0]['n_stars'],
      scenario_mean=sc.round(4).to_dict(),top_ID_scenario=[(int(i),s,round(float(v),4)) for (i,s),v in top.items()])
json.dump(out,open('evidence/05_triceratops_summary.json','w'),indent=1)
for d,v in out.items(): print(d,json.dumps(v))
