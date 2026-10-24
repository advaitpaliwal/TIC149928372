# TIC v8.2 query for hz-bright: 0<d<100 pc, 2700<=Teff<=4000, T<=13.5, STAR (logg/R/M cuts applied in build_pool.py so blank logg survives)
from pathlib import Path
from astroquery.mast import Catalogs
r=Path(__file__).resolve().parent;out=[]
import astropy.table as at
for lo,hi in [(-5,9),(9,10.5),(10.5,11.5),(11.5,12.2),(12.2,12.8),(12.8,13.2),(13.2,13.5)]:
    q=Catalogs.query_criteria(catalog='Tic',d=[0,100],Teff=[2700,4000],Tmag=[lo,hi],objType='STAR')
    q=q[['ID','ra','dec','Tmag','Teff','rad','mass','d','logg','rho','lumclass','disposition','wdflag']];print(lo,hi,len(q),flush=True);out.append(q)
t=at.vstack(out);t.write(r/'_sources/tic/tic_d100.csv',format='ascii.csv',overwrite=True);print(len(t))
