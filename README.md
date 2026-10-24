# TIC 149928372.01: a 2.7 Earth-radius, 16.36-day transiting planet candidate around an early M dwarf

Discovered by Advait Paliwal on 2026-10-23 in public TESS photometry (search result written 2026-10-23 21:47 PDT, `data/search/149928372_search_result.json`).

## Summary

TIC 149928372 (Gaia DR3 1058588908797940864) is an early M dwarf (Teff 4000 ± 125 K, 0.66 R☉, 0.62 M☉) at 79 pc. TESS light curves from Sectors 14 and 21 (SPOC 2-minute) and 41, 47, 48, 74 and 75 (TESS-SPOC full-frame images), spanning 2019 to 2024, show 13 transits with period 16.36448 d, depth about 1,450 ppm and duration 3.36 h, which implies a radius of 2.67 (+0.29 / −0.28) Earth radii and an insolation of 8.7 times Earth's. The transit is flat-bottomed (impact parameter 0.37, no grazing posterior samples), and the stellar density fitted from the transit shape, 3.33 (+0.74 / −1.42) g cm⁻³, matches the catalogue value of 3.04 ± 1.14. Odd and even depths agree (1.25σ), there is no secondary eclipse, the dip is in the raw (SAP) flux (1,443 ± 68 against 1,450 ± 74 ppm), and difference images in all seven sectors place it on the target. As of 2026-10-24 the star is not a TOI or a CTOI, and it is not a NASA Exoplanet Archive host.

**This is a planet candidate, not a validated or confirmed planet.** With no neighbours removed, the TRICERATOPS false-positive probability is 0.0169 ± 0.0035 and the nearby false-positive probability is 1.5 × 10⁻⁹. With Gaia companion limits from MOLUSC (Wood, Mann & Kraus 2021, AJ 162, 128) it is 0.0019 ± 0.0017. High-resolution imaging has not been checked (ExoFOP was unreachable on 2026-10-24), so no validation claim is made.

The signal was a SPOC full-frame-image threshold-crossing event in Sector 47 (`149928372-1-S47`, P 16.3779 d, MES 8.8, two transits). Its data-validation report put the difference-image centroid 76.1 ± 1.0″ from the target, toward the G = 8.9 giant TIC 149928367 at 83″, and ExoMiner++ 2.0 (Valizadegan et al. 2025; Martinho & Valizadegan 2025, Zenodo 17707413) scored it 0.037, a false positive. Our seven-sector pixel fits put the dip on the target: the fitted flux loss at the target is 0.87 ± 0.13 of what the light-curve depth requires, and forcing the signal onto the giant is excluded at 8 to 63σ per sector. In our Sector 47 difference image the giant's fitted amplitude is 2.4σ against 97 random-epoch null images, so the Sector 47 offset is not significant.

| | |
|---|---|
| Star | TIC 149928372, RA 162.06691°, Dec +65.10011° (J2000), Tmag 11.40, G 12.22 |
| Period | 16.364477 ± 0.000035 d |
| Mid-transit (T0) | BJD_TDB 2459518.8453 ± 0.0013 |
| Depth | 1,450 ± 74 ppm |
| Duration | 3.36 h |
| Impact parameter | 0.37 (+0.27 / −0.25) |
| Radius | 2.67 (+0.29 / −0.28) R⊕ |
| Insolation | 8.7 S⊕ |
| False-positive probability | 0.0169 ± 0.0035 (no companion limits), 0.0019 ± 0.0017 (MOLUSC Gaia limits); no neighbours removed |
| High-resolution imaging | not checked |

## Contents

- `note/TIC149928372_note.md`, `note/TIC149928372_note.pdf`: the research note (draft; not submitted to a journal).
- `note/figure1_TIC149928372.png`: the seven sectors folded on the period, the 13 transit depths, odd and even transits, and phase 0.5.
- `ephemeris.csv`: predicted mid-transit times from 2026-10-24 to 2027-12-31 with 1σ timing uncertainty.
- `code/pool/`: how the 31,405-star search pool was built (`pool_counts.json` has the count at each cut).
- `code/search/`: the star search (`run1.py` light-curve reader for this pool, `run_star_gpu.py` and `drive_gpu.py` for the cuvarbase GPU box search, `GPU.md` for its calibration and control planets, `run_star.py` for the earlier Transit Least Squares search, `edgefix/` for the light-curve preparation tests).
- `code/leo/`: the LEO-vetter gate (`gate.py`, `leo_glue.py`, `control.py`, `LEO.md`).
- `code/vet/`: the vetting scripts (native-cadence windows, transit fit, shape and secondary checks, SAP against PDCSAP, pixel difference images and null tests, TRICERATOPS and MOLUSC runs).
- `code/ephemeris/`: the posterior conversion (`make_chain.py`) and the prediction table (`ephemeris.py`).
- `code/note/`: Figure 1 (`fig1.py`) and the Word build of the note.
- `data/search/`: the search result for this star and the LEO-vetter control results.
- `data/triage/`: triage of the first 14 stars from this search to pass LEO-vetter (this one, 4 known signals, 9 false positives).
- `data/vet/`: the vetting report (`VET.md`), its evidence files and plots, the flattened search light curve (`149928372.npz`), the saved SPOC Sector 47 data-validation report and the MAST URLs of the light curves.
- `data/pixels/`: stacked difference images per sector, the Gaia field and the neighbour tests.
- `data/triceratops/`, `data/triceratops_molusc/`: per-run TRICERATOPS outputs without and with MOLUSC limits.
- `data/molusc/`: MOLUSC logs and summary.
- `data/fit/chain.npz`: the ephemeris posterior.
- `data/exofop_ctoi_row.txt`: the prepared ExoFOP community-candidate parameter row (not uploaded).

## Next predicted transits

From `ephemeris.csv` (1σ timing uncertainty 5.6 minutes now, growing to 6.8 minutes by end of 2027):

| Epoch | Mid-transit BJD_TDB | Mid-transit UTC | ±1σ |
|---|---|---|---|
| 112 | 2461351.66668 | 2026-11-07 03:58 | 5.6 min |
| 113 | 2461368.03116 | 2026-11-23 12:41 | 5.6 min |
| 114 | 2461384.39564 | 2026-12-09 21:25 | 5.7 min |
| 115 | 2461400.76012 | 2026-12-26 06:09 | 5.7 min |
| 116 | 2461417.12460 | 2027-01-11 14:53 | 5.7 min |
| 117 | 2461433.48907 | 2027-01-27 23:38 | 5.8 min |

Times are propagated from every posterior sample of the vetting transit fit, in which the period and mid-transit time are both free (`code/vet/s1_fit.py free`), so the correlation between T0 and P is kept. UTC times are geocentric (BJD_TDB minus the light-travel delay to the star, then TDB to UTC); ingress and egress are mid-transit ± 1.68 h. Epochs count from T0 = BJD_TDB 2459518.8453.

## Reproducing

Python 3.12 with numpy, scipy, pandas, matplotlib, astropy, astroquery, lightkurve, wotan, cuvarbase (GPU search), batman-package, emcee, TESS_PRF, LEO-vetter 1.2.0 (https://github.com/mkunimoto/LEO-Vetter), TRICERATOPS 1.1.0 and MOLUSC (https://github.com/woodml/MOLUSC, not bundled).

1. Light curves. SPOC 2-minute PDCSAP for Sectors 14 and 21 and TESS-SPOC PDCSAP for Sectors 41, 47, 48, 74 and 75, from MAST (URLs in `data/vet/lcfits_urls/`). FITS files and TESScut cutouts are not included; they are public at MAST. The flattened search light curve is `data/vet/149928372.npz`; the native-cadence transit windows are `data/vet/evidence/01a_windows.npy`.
2. Search. `code/search/run_star_gpu.py` (reads light curves through `run1.py`), then `code/leo/gate.py` for LEO-vetter.
3. Vetting. `code/vet/` from a working folder holding `evidence/`: `s1_prep.py`, `s1_fit.py free`, `s1_fit.py prior`, `s1_plots.py`, `s1_shape.py`, `sap_pdc.py`, `pix_dl.py`, `diff_all.sh`, `post.py`, `pix_null.py`, then `trc_run.py` and `trc_summary.py`, then `molusc_all.sh`, `molusc_summary.py` and `trc_run.py` with the MOLUSC file. The 10⁶-row MOLUSC file passed to TRICERATOPS and the full MOLUSC outputs are not included because of their size; `molusc_summary.py` regenerates them from MOLUSC.
4. Ephemeris. `python code/ephemeris/make_chain.py data/vet/evidence/01b_chain_free.npy data/fit/chain.npz`, then `python code/ephemeris/ephemeris.py data/fit/chain.npz ephemeris.csv`.

Some scripts use paths relative to the original working folders; adjust them to your layout.

## Credits

This work uses data collected by the TESS mission, funded by the NASA Science Mission Directorate, obtained from the Mikulski Archive for Space Telescopes (MAST), including SPOC (Jenkins et al. 2016) and TESS-SPOC (Caldwell et al. 2020) light curves, the SPOC data-validation report and TESScut cutouts. It uses data from the European Space Agency mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium, and the TESS Input Catalog (Stassun et al. 2019). Companion limits use MOLUSC (Wood, Mann & Kraus 2021, AJ 162, 128); false-positive probabilities use TRICERATOPS (Giacalone & Dressing 2020; Giacalone et al. 2021); vetting uses LEO-vetter (Kunimoto et al. 2025). Full references are in the note.

## Licence and citation

Code is under the MIT licence (`LICENSE`). Text, figures and data are under CC BY 4.0 (`LICENSE-CC-BY-4.0`). Please cite with `CITATION.cff`.

## Contact

Advait Paliwal, advaitpaliwal@gmail.com
