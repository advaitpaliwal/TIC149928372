# TOI-700 (TIC 150428135) control: TIC row + MAST products, through the same build_pool.py product rule -> control.csv
from pathlib import Path
import pandas as pd
from astroquery.mast import Catalogs,Observations
r=Path(__file__).resolve().parent
q=Catalogs.query_criteria(catalog='Tic',ID=150428135)[['ID','ra','dec','Tmag','Teff','rad','mass','d','logg','rho']].to_pandas()
q.to_csv(r/'control_pre.csv',index=False)
o=Observations.query_criteria(target_name='150428135',obs_collection=['TESS','HLSP'],dataproduct_type='timeseries',provenance_name=['SPOC','TESS-SPOC','QLP'])
o[['target_name','provenance_name','sequence_number','t_exptime','dataURL']].to_pandas().to_csv(r/'control_obs.csv',index=False);print(q)
