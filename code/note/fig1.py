# Figure 1 for the TIC 149928372 note.
# Usage: python -I fig1.py <01a_windows.npy> <01b_fit_free.json> <01c_epochs_oddeven.json> <149928372.npz> <01d_shape_secondary.json> <out.png>
import sys, json, numpy as np, batman
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

A = np.load(sys.argv[1]); t, f, e, ex, sec, ep = A.T          # native-cadence PDCSAP windows (vet glue/s1_prep.py)
fit = json.load(open(sys.argv[2])); oe = json.load(open(sys.argv[3]))
d = np.load(sys.argv[4]); ts, fs = d['t'], d['f']              # flattened search light curve, 10-min
sec2 = json.load(open(sys.argv[5]))['secondary_phase05']
P, T0 = fit['P'][1], fit['T0'][1]                              # free-fit posterior medians
rp, b, a = fit['rp_rs'][1], fit['b'][1], fit['a_rs'][1]

def model(x):
    pr = batman.TransitParams(); pr.t0 = 0; pr.per = P; pr.rp = rp; pr.a = a
    pr.inc = np.degrees(np.arccos(b/a)); pr.ecc = 0; pr.w = 90; pr.limb_dark = 'quadratic'; pr.u = fit['LD_fixed']
    return batman.TransitModel(pr, x, supersample_factor=11, exp_time=10/1440).light_curve(pr)

def binned(x, y, w):
    edges = np.arange(x.min(), x.max()+w, w); i = np.digitize(x, edges)
    k = [k for k in range(1, len(edges)) if (i == k).sum() > 2]
    return (np.array([(edges[j-1]+edges[j])/2 for j in k]), np.array([np.mean(y[i == j]) for j in k]),
            np.array([np.std(y[i == j])/np.sqrt((i == j).sum()) for j in k]))

dt = (((t-T0)/P+0.5) % 1 - 0.5)*P; n = np.round((t-T0)/P).astype(int)
ppt = lambda y: (y-1)*1e3
ink, grey, accent = 'k', '0.75', '#c4552b'
xm = np.linspace(-0.3, 0.3, 1200)
fig = plt.figure(figsize=(7.5, 8.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 0.9, 0.9], hspace=0.55, wspace=0.3)

ax = fig.add_subplot(gs[0, :]); w = np.abs(dt) < 0.3
ax.plot(dt[w]*24, ppt(f[w]), '.', color=grey, ms=1.2, zorder=1)
c, m, r = binned(dt[w], f[w], 10/1440)
ax.errorbar(c*24, ppt(m), r*1e3, fmt='o', color=ink, ms=3, lw=0.8, zorder=3)
ax.plot(xm*24, ppt(model(xm)), '-', color=accent, lw=1.4, zorder=4)
ax.set_xlim(-7, 7); ax.set_ylim(-3, 1.5); ax.set_xlabel('Hours from mid-transit'); ax.set_ylabel('Relative flux (ppt)')
ax.set_title(f'(a) 7 sectors folded on P = {P:.5f} d, 10-min bins', fontsize=9, loc='left')

ax = fig.add_subplot(gs[1, :]); pe = oe['per_epoch']
x = np.array([p['epoch'] for p in pe]); y = np.array([p['depth_ppm'] for p in pe]); s = np.array([p['err_ppm'] for p in pe])
ax.errorbar(x, y, s, fmt='o', color=ink, ms=3.5, lw=0.8)
ax.axhline(oe['weighted_mean_ppm'], color=accent, lw=1, ls='--')
seen = set()
for p in pe:
    if p['sector'] in seen: continue
    seen.add(p['sector']); ax.annotate(f"S{p['sector']}", (p['epoch'], p['depth_ppm']+p['err_ppm']+60), fontsize=6.5, ha='center')
ax.set_ylim(0, 3000); ax.set_xlabel('Transit number'); ax.set_ylabel('Depth (ppm)')
ax.set_title(f"(b) the 13 transits; $\\chi^2$ = {oe['chi2']:.1f} for {oe['dof']} dof about the mean (dashed)", fontsize=9, loc='left')

ax = fig.add_subplot(gs[2, 0])
for sel, col, lab, off in [(n % 2 == 1, ink, 'odd', 0), (n % 2 == 0, accent, 'even', 0.08)]:
    c, m, r = binned(dt[w & sel], f[w & sel], 20/1440)
    ax.errorbar(c*24+off, ppt(m), r*1e3, fmt='o', color=col, ms=2.5, lw=0.7, label=lab)
ax.set_xlim(-7, 7); ax.set_ylim(-3, 1.5); ax.legend(fontsize=7, frameon=False, loc='lower right')
ax.set_xlabel('Hours from mid-transit'); ax.set_ylabel('Rel. flux (ppt)')
ax.set_title('(c) odd and even transits', fontsize=9, loc='left')

ax = fig.add_subplot(gs[2, 1]); d5 = (((ts-T0-P/2)/P+0.5) % 1 - 0.5)*P; w5 = np.abs(d5) < 0.3
ax.plot(d5[w5]*24, ppt(fs[w5]), '.', color=grey, ms=1.5)
c, m, r = binned(d5[w5], fs[w5], 30/1440)
ax.errorbar(c*24, ppt(m), r*1e3, fmt='o', color=ink, ms=2.5, lw=0.7)
ax.set_xlim(-7, 7); ax.set_ylim(-3, 1.5); ax.set_xlabel('Hours from phase 0.5')
ax.set_title('(d) phase 0.5, 30-min bins', fontsize=9, loc='left')

fig.savefig(sys.argv[6], dpi=200, bbox_inches='tight')
print('model depth ppm', (1-model(np.array([0.0]))[0])*1e6, 'n windows pts', w.sum(), 'phase0.5 pts', w5.sum(), 'sec0.5', sec2)
