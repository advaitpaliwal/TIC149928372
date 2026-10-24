# usage: build_pool.py [IN OBS OUT]  (default pre_products.csv _sources/mast/obs.csv pool.csv)
# hz-near's product rule: one light curve per sector, SPOC 120-s first, then TESS-SPOC, then QLP; keep >=3 sectors.
# Output sorted brightest first. Appends counts to pool_counts.json when building the pool.
from pathlib import Path
import pandas as pd,numpy as np,json,sys
r=Path(__file__).resolve().parent
a=sys.argv[1:] or ['pre_products.csv','_sources/mast/obs.csv','pool.csv']
t=pd.read_csv(r/a[0]);o=pd.read_csv(r/a[1]);o['tic']=o.target_name.astype(int)
o=o[o.dataURL.astype(str).str.endswith('lc.fits')]
o=o[((o.provenance_name=='SPOC')&(o.t_exptime.round()==120))|(o.provenance_name.isin(['TESS-SPOC','QLP']))]
o['rk']=o.provenance_name.map({'SPOC':0,'TESS-SPOC':1,'QLP':2})
o=o.sort_values(['tic','sequence_number','rk']).drop_duplicates(['tic','sequence_number'])
def summ(d):return pd.Series({'sectors':json.dumps([[int(s),p,float(e),u] for s,p,e,u in zip(d.sequence_number,d.provenance_name,d.t_exptime,d.dataURL)]),'ns':len(d),'n_2min':int((d.provenance_name=='SPOC').sum()),'n_tess_spoc':int((d.provenance_name=='TESS-SPOC').sum()),'n_qlp':int((d.provenance_name=='QLP').sum())})
g=o.groupby('tic').apply(summ).reset_index().rename(columns={'tic':'ID'})
t=t.merge(g,on='ID',how='left').fillna({'ns':0,'n_2min':0,'n_tess_spoc':0,'n_qlp':0})
st={'into_product_cut':len(t),'with_any_product':int((t.ns>0).sum()),'with_1_or_2_sectors':int(t.ns.between(1,2).sum()),'with_0':int((t.ns==0).sum())}
t=t[t.ns>=3].sort_values('Tmag');t['ns']=t.ns.astype(int)
st.update(pool=len(t),pool_with_any_2min=int((t.n_2min>0).sum()),pool_2min_only=int((t.n_2min==t.ns).sum()),pool_Tmag_lt_11_5=int((t.Tmag<11.5).sum()),pool_Tmag_11_5_to_13_5=int((t.Tmag>=11.5).sum()),
        pool_teff_outside_hunt7_2800_3800=int((~t.Teff.between(2800,3800)).sum()),pool_blank_logg=int(t.logg.isna().sum()),
        products_2min=int(t.n_2min.sum()),products_tess_spoc=int(t.n_tess_spoc.sum()),products_qlp=int(t.n_qlp.sum()),median_sectors=float(t.ns.median()))
t.drop(columns=['lumclass','disposition','wdflag'],errors='ignore').to_csv(r/a[2],index=False);print(json.dumps(st,indent=1))
if a[2]=='pool.csv':
    c=json.load(open(r/'counts1.json'));c.update(st);json.dump(c,open(r/'pool_counts.json','w'),indent=1)
