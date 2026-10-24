# hz-bright cuts before product lookup. Writes pre_products.csv and counts1.json.
from pathlib import Path
import pandas as pd,numpy as np,json
from astropy.coordinates import SkyCoord
import astropy.units as u
r=Path(__file__).resolve().parent;S=r/'_sources'
t=pd.read_csv(S/'tic/tic_d100.csv');st={'tic_query_rows':len(t)}
t=t.drop_duplicates('ID').copy();t.ID=t.ID.astype(int);st['unique_ids']=len(t)
t=t[t.d.gt(0)&t.d.lt(100)&t.Teff.between(2700,4000)&t.Tmag.le(13.5)];st['d_teff_tmag']=len(t)
t=t[t.disposition.isna()];st['no_tic_disposition(dup/split/artifact)']=len(t)
t=t[t.lumclass.ne('GIANT')];st['not_giant']=len(t)
t=t[t.logg.isna()|t.logg.gt(4.3)];st['logg_gt_4.3_or_blank']=len(t);st['  of_which_blank_logg']=int(t.logg.isna().sum())
t=t[t.rad.gt(0)&t.mass.gt(0)];st['R_and_M_positive']=len(t)
sr=pd.read_csv(S/'searched_ids.csv',header=None,names=['lane','dir','ID']);done=set(sr.ID.astype(int))
for l in ['hz-north','hz-south','hz-near','hz-wide']:st[f'  in_{l}']=int(t.ID.isin(set(sr[sr.lane==l].ID)).sum())
t=t[~t.ID.isin(done)];st['minus_searched_any_lane']=len(t)
toi=set(pd.read_csv(S/'cat/toi.csv')['TIC ID'].dropna().astype(int))
ctoi=set(pd.read_csv(S/'cat/ctoi.csv')['TIC ID'].dropna().astype(int))
nea=set(pd.read_csv(S/'cat/ps_tic.csv').tic_id.dropna().str.replace('TIC ','').astype(int))
ps=pd.read_csv(S/'cat/ps_coords.csv').dropna(subset=['ra','dec'])
c=SkyCoord(t.ra.values*u.deg,t.dec.values*u.deg);pc=SkyCoord(ps.ra.values*u.deg,ps.dec.values*u.deg)
_,sep,_=c.match_to_catalog_sky(pc);near=sep.arcsec<20
f={'toi':t.ID.isin(toi).values,'ctoi':t.ID.isin(ctoi).values,'nea_tic':t.ID.isin(nea).values,'nea_within_20arcsec':near}
for k,v in f.items():st['  excl_'+k]=int(v.sum())
bad=np.any(np.vstack(list(f.values())),axis=0);t=t[~bad];st['minus_toi_ctoi_nea_hosts']=len(t)
st['  note']='SPOC TCE hosts are NOT dropped'
t=t.sort_values('Tmag');t.to_csv(r/'pre_products.csv',index=False);json.dump(st,open(r/'counts1.json','w'),indent=1);print(json.dumps(st,indent=1))
