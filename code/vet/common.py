import warnings; warnings.filterwarnings('ignore')
import numpy as np, glob, json, lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
_c=json.load(open('evidence/01b_fit_free.json'))
P,T0=_c['P'][1],_c['T0'][1]; DUR=_c['T14_h'][1]/24
RA0,DE0=162.066908906837,65.1001060259593; GT=12.2324
def tpf_for(s):
    p=glob.glob(f'pix/cut/tesscut/tess-s{s:04d}-*_21x21_astrocut.fits')[0]
    return lk.read(p)
def gaia_field(tpf, rad=200):
    import os, pandas as pd
    fn='pix/gaia_field.csv'
    if not os.path.exists(fn):
        v=Vizier(columns=['Source','RA_ICRS','DE_ICRS','pmRA','pmDE','Gmag','BPmag','RPmag','Plx','RUWE'],row_limit=-1)
        t=v.query_region(SkyCoord(RA0*u.deg,DE0*u.deg),radius=rad*u.arcsec,catalog='I/355/gaiadr3')[0].to_pandas()
        t.to_csv(fn,index=False)
    t=pd.read_csv(fn)
    dy=(tpf.time[0].jd-2457389.0)/365.25
    ra=(t.RA_ICRS+np.nan_to_num(t.pmRA)*dy/3.6e6/np.cos(np.radians(t.DE_ICRS))).values; de=(t.DE_ICRS+np.nan_to_num(t.pmDE)*dy/3.6e6).values
    x,y=tpf.wcs.all_world2pix(np.vstack([ra,de]).T,0).T
    c=SkyCoord(ra*u.deg,de*u.deg); i0=int(np.argmin(np.nan_to_num(np.abs(t.Gmag.values-GT),nan=99.)+(c.separation(SkyCoord(RA0*u.deg,DE0*u.deg)).arcsec>20)*99))
    t=t.assign(ra_ep=ra,de_ep=de,x=x,y=y,sep=c.separation(c[i0]).arcsec,pa=c[i0].position_angle(c).deg)
    return t.sort_values('sep').reset_index(drop=True)
