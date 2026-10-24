# MAST light-curve products (SPOC 2-min, TESS-SPOC, QLP) for every star in pre_products.csv, threaded. Writes _sources/mast/obs.csv
from pathlib import Path
import pandas as pd,sys,time
from concurrent.futures import ThreadPoolExecutor
from astroquery.mast import Observations
r=Path(__file__).resolve().parent
ids=[str(int(i)) for i in pd.read_csv(r/'pre_products.csv').ID];B=300
def q(i):
    b=ids[i:i+B]
    for k in range(5):
        try:
            o=Observations.query_criteria(target_name=b,obs_collection=['TESS','HLSP'],dataproduct_type='timeseries',provenance_name=['SPOC','TESS-SPOC','QLP'])
            print(i,len(o),flush=True)
            return o[['target_name','provenance_name','sequence_number','t_exptime','dataURL']].to_pandas() if len(o) else None
        except Exception as e:print('retry',i,repr(e)[:120],flush=True);time.sleep(10*(k+1))
    raise SystemExit(f'batch {i} failed')
with ThreadPoolExecutor(8) as ex:out=[d for d in ex.map(q,range(0,len(ids),B)) if d is not None]
pd.concat(out).to_csv(r/'_sources/mast/obs.csv',index=False);print('done',sum(map(len,out)))
