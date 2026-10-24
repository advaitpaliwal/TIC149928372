import sys, lightkurve as lk, warnings; warnings.filterwarnings('ignore')
for s in map(int,sys.argv[1:]):
    try:
        r=lk.search_tesscut('TIC 149928372',sector=s); tpf=r.download(cutout_size=21,download_dir='pix/cut')
        print(s,tpf.path, tpf.shape, flush=True)
    except Exception as z: print(s,'ERR',z,flush=True)
