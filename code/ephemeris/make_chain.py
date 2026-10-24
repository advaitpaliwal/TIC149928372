# Posterior for the ephemeris of TIC 149928372.01. The vetting fit (code/vet/s1_fit.py free) already
# has T0 and P free, so its chain carries the T0-P correlation; this script only converts it.
# Usage: python make_chain.py 01b_chain_free.npy chain.npz   (columns T0 [BTJD], P, Rp/R*, b, ln a/R*)
import sys, numpy as np
ch = np.load(sys.argv[1])
T0 = ch[:, 0] + 2457000.0; P = ch[:, 1]
np.savez(sys.argv[2], T0_BJD=T0, P=P, rp_rs=ch[:, 2], b=ch[:, 3], a_rs=np.exp(ch[:, 4]))
q = lambda v: np.percentile(v, [16, 50, 84])
print('n', len(P), 'T0', q(T0), 'std', T0.std(), 'P', q(P), 'std', P.std(), 'corr', np.corrcoef(T0, P)[0, 1])
