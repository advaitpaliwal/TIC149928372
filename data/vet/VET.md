# TIC 149928372 vet (hz-bright GPU search)

Run 2026-10-24, 02:12–02:42 EDT. The Mac mini work (`hunt4/.venv`, `python -I`, `taskpolicy -a`, at most 3 single-thread jobs) covered the light-curve checks, pixels and MOLUSC. The transit MCMC and TRICERATOPS ran on the Spark, in `~/planet-hunt/trc_149928372` with `trc_venv`, `nice 19`, `NUMPY_MADVISE_HUGEPAGE=0` and `NUMBA_NUM_THREADS=1`. I moved them there because the Mac mini load average was 30–40. The GPU search and `~/planet-hunt/.venv` were not touched. All FITS files were deleted after reading, including the Spark lightkurve TESScut cache files for this star. Their URLs are kept as `.url` files.

Star (TIC 8.2, `_sources/tic82_149928372.ecsv`): K/M dwarf, Teff 4000 ± 125 K, R 0.66 ± 0.07 R☉, M 0.62 ± 0.08 M☉, density 2.16 ± 0.81 solar (3.04 ± 1.14 g/cc), Tmag 11.40, 79.4 pc. Gaia DR3 1058588908797940864, G 12.22, RUWE 1.16 (`pix/gaia_field.csv`).

## Verdict

**CANDIDATE.**

- **TRICERATOPS with no neighbours removed, 10 runs: FPP 0.0169 ± 0.0035 and NFPP 1.5e-9.** This passes the candidate bar (FPP ≤ 0.035) and the NFPP bar (≤ 0.001). It misses the validated-grade bar (≤ 0.015) by a little: 4 of the 10 runs are under 0.015 (`evidence/05_triceratops_summary.json`).
- **With MOLUSC Gaia companion limits, 10 runs: FPP 0.0019 ± 0.0017 and NFPP 1.1e-9.** That is validated-grade.
- The grazing / V-shape worry from triage does not hold up on the native-cadence data:
  - The dip is flat-bottomed.
  - The impact parameter is b ≈ 0.37.
  - The density fitted from the transit shape matches the TIC star.
- The planet comes out bigger than triage said: **Rp ≈ 2.65 R⊕, depth about 1,450 ppm.** Triage had 2.2 R⊕ and 903 ppm.
- What is left in the FPP is "planet on a bound companion": PTP 0.14 and STP 0.017 of the probability. With MOLUSC, STP falls to 0.002. Without high-resolution imaging, nothing removes the PTP term.

## 1. Transit fit: the crux

Method:
- Native-cadence PDCSAP from all 7 sectors: SPOC 2-min for S14 and S21, TESS-SPOC 10-min for S41, S47 and S48, and 200-s for S74 and S75.
- I cut a ±0.6 d window around each of the 13 transits and divided out a quadratic fit outside ±0.15 d (`glue/s1_prep.py`, `evidence/01a_windows.json`).
- Fit: batman + emcee with rp, b, a/R★, T0 and P free. Limb darkening was fixed at u1 0.44 and u2 0.20, my approximate value for 4000 K (not looked up). Each cadence is supersampled (`glue/s1_fit.py`).

| fit | b | a/R★ | ρ★ (g/cc) | Rp/R★ | Rp (R⊕) | T14 (h) | T23/T14 | P(grazing) |
|---|---|---|---|---|---|---|---|---|
| free density, all data | 0.37 (+0.27/−0.25) | 36.1 (+2.5/−6.1) | 3.33 (+0.74/−1.42) | 0.0369 (+0.0017/−0.0012) | 2.66 (2.38–2.96 with the R★ error) | 3.36 | 0.92 | 0.000 |
| TIC density prior | 0.36 (+0.20/−0.23) | 36.2 | 3.36 (+0.64/−0.98) | 0.0368 | 2.65 | 3.36 | 0.92 | 0.000 |
| TIC expectation | – | 35.0 | 3.04 ± 1.14 | – | – | – | – | – |

[verified, `evidence/01b_fit_free.json`, `evidence/01b_fit_prior.json`, `plots/01_transit_fit.png`]

- With the density left free, the fit returns 3.33 g/cc. That is 0.3σ from the TIC density.
- Only 2.8% of the posterior has b > 0.8, and none is grazing. Only 1.5% has Rp > 3 R⊕.
- Ephemeris: **P 16.364477 ± 0.000033 d, T0 BTJD 2518.8453 ± 0.0012.**
- The 2-min data alone (S14+S21) give the same U shape (right panel of `plots/01_transit_fit.png`).

**Forcing a grazing fit** [verified, `evidence/01d_shape_secondary.json`]. b was fixed and the other parameters refit. Δχ² is measured against the free fit (χ² 6,210 on 5,930 points).

| b fixed | Δχ² | Rp (R⊕) | ρ★ implied (g/cc) |
|---|---|---|---|
| 0.85 | +6.3 | 3.0 | 0.66 |
| 0.90 | +10.1 | 3.2 | 0.39 |
| 0.95 | +20.1 | 3.5 | 0.16 |
| 1.00 | +24.7 | 5.3 | 0.09 |

A grazing solution fits worse. It would also need a star 8 to 30 times less dense than a 0.62 M☉ dwarf. [inferred] Triage's LEO fit (b 0.90, a/R★ 16) and its "between U and V" shape came from the 30-minute, flattened search light curve. That light curve also gave a shallower depth (903 ppm versus about 1,450 ppm here), so the search's detrending probably ate part of the dip.

**Every transit and odd/even** [verified, `evidence/01c_epochs_oddeven.json`, `plots/01_per_epoch.png`, `plots/01_every_epoch.png`]:
- 13 transits, weighted mean 1,469 ppm, χ² 19.2 on 12 degrees of freedom (p ≈ 0.08).
- Odd 1,381 ± 103 ppm and even 1,572 ± 112 ppm, 1.25σ apart.

**Secondary eclipse** [verified, `evidence/01d_shape_secondary.json`, `plots/01_secondary_phases.png`]:
- Phase 0.5 measures −133 ± 114 ppm, so no secondary eclipse.
- The largest value on a 91-phase grid is 421 ± 113 ppm (3.7σ) at phase 0.54. The grid's scatter is 1.7σ, so this is about 2.2 times the scatter, which is noise-like [inferred].

**SPOC's own S47 DV fit** [verified, `evidence/02_spoc_s47_dv.json`, `_sources/*dvr.xml`]:
- b 0.02, a/R★ 33.1, stellar density 1.82 solar, Rp 2.97 ± 1.44 R⊕.
- `suspectedEclipsingBinary="false"`, odd/even significance 0.75, weak secondary 580 ± 243 ppm.

## 2. Pixels

**Method.** TESScut 21×21 cutouts for all 7 sectors. For each sector I built an in-transit minus out-of-transit difference image and calibrated it against 120 random-epoch null images. I fit the TESS PRF at the target and at every Gaia star within 150″ down to G 20.5. The glue was copied from `hz-wide/vet_270573718`: `glue/common.py` differs only in the ephemeris source and coordinates, and `diff.py` is unchanged.

**The target matches** [verified, `evidence/03_target_amp.csv`, `evidence/03b_pixel_null_joint.json`]:
- In every sector, the PRF amplitude at the target matches the amplitude the light-curve depth requires.
- The weighted ratio of fitted to required amplitude is **0.87 ± 0.13**.
- S75 is the low one: 0.38, 1.8σ low. It has only one transit, next to a data gap.

**Neighbours within 2′.** Only two are bright enough to matter. Fainter stars, out to 150″, would need a deeper-than-100% eclipse.
- **G 19.1 at 29″ (TIC 149928373):** it would need a 61–113% eclipse, depending on the sector. A forced fit is rejected in every sector: summed Δχ² 1,720, and null-calibrated per-sector exclusions of 2.7 to 15.5σ.
- **G 8.9 giant at 83″ (TIC 149928367, Tmag 8.2, R 11 R☉):** it would need a 0.4–3.5% eclipse. A forced fit is rejected in every sector, at 8.0 to 63σ calibrated against the nulls.
- The other stars within 2′: G 21.0 at 34″, G 20.3 at 53″, G 20.2 at 68″ and G 21.8 at 71″. Each would need more than 100% (`evidence/03_neighbour_summary.csv`).
  - Note: that file's "ruled out by pixels" label on the target row (Δχ² 33) comes from the low S21 and S75 amplitudes. It is not a neighbour result.

**Why ExoMiner++ scored it a false positive** [verified data, inferred cause]:
- SPOC's S47 difference-image centroid sits **76.1 ± 1.0″ from the target**, toward the giant (`ticCentroidOffsets` in the DV xml).
- My S47 difference image also shows a dimming blob at the giant (`plots/03_diff_s47.png`). Its free amplitude is 100 e/s.
- Calibrated against the null images, though, that blob is only **2.4σ**. The giant's amplitudes in the other sectors are −1.1 to +1.0σ.
- The giant would have needed about 600 e/s in S47 to make the dip (11.8σ excluded).
- [inferred] The red giant's own brightness wiggles leaked into a 2-transit difference image. That pulled SPOC's centroid toward it, and this is the likely source of the 0.037 score.

## 3. SAP vs PDCSAP

[verified, `evidence/04_sap_vs_pdcsap.json`, `plots/04_sap_vs_pdcsap.png`; local linear baseline per event, no detrending]
- The dip is in the raw flux.
- All 7 sectors: SAP 1,443 ± 68 ppm, SAP/CROWDSAP 1,623 ± 77 ppm, PDCSAP 1,450 ± 74 ppm.
- Without S75: SAP/CROWDSAP 1,460 ± 79 ppm against PDCSAP 1,432 ± 78 ppm.
  - S75's single raw event (4,828 ppm) sits on a data-gap ramp. Its PDCSAP value is 1,633 ppm.
- The depths agree across sectors even though CROWDSAP ranges from 0.75 (S21) to 0.96 (S14). [inferred] That is what you expect if the source is the target, not a contaminant.

## 4. TRICERATOPS

Setup [verified, `glue/trc_run.py`, `trc/out`, `trc_molusc/out`]:
- 7 sectors, the SPOC/TESS-SPOC pipeline apertures, and the native-cadence windows folded into 10-min bins over ±0.4 d. Depth 1,428 ppm, bin error 417 ppm.
- 33 stars in the field. Three have a non-zero required depth: the target, the G 19.1 star (0.70) and the giant (0.24).
- N = 10⁶ samples, seeds 1000–1009.

| run | FPP | NFPP | where the false-positive probability sits |
|---|---|---|---|
| no neighbours removed | **0.0169 ± 0.0035** (0.0100–0.0207) | 1.5e-9 | PTP 0.142, STP 0.017, DTP 0.001; every NEB/NTP is 0 |
| + MOLUSC companion limits | **0.0019 ± 0.0017** (0–0.0054) | 1.1e-9 | PTP 0.137, STP 0.002 |

MOLUSC [verified, `molusc_out/summary.txt`]:
- 4 runs of 50,000 companions, using the RUWE and Gaia contrast limits.
- 50.4% of companions survive. RUWE rejects 19.5% and Gaia contrast rejects 31.2%.
- The file passed to TRICERATOPS is `molusc_out/molusc_for_trc_N1e6.csv`.

No neighbour removal was needed, because the nearby-star terms are already zero.

## 5. Catalogues and ExoFOP

- **ExoFOP:** still down. Direct curl from the Mac mini and the Spark gave Cloudflare 525 and then timeouts, and Exa's fetch timed out too (`_sources/exofop_attempts.txt`).
  - *Not checked:* imaging, spectroscopy and observing notes.
  - Triage found no TOI or CTOI (`../triage1/TRIAGE.md`).
- **Only one SPOC product exists, the S47 FFI DV report:** P 16.3779 d, MES 8.85, 2 transits (`_sources/*dvr.xml`, `*dvs-01.pdf`).

## Limits

- There is no high-resolution imaging and no spectroscopy. The PTP term (a planet on an unresolved bound companion) is about 0.14 of the probability, and only imaging can remove it.
- Limb darkening was fixed at approximate values, not taken from Claret tables.
- The orbit was assumed circular.
- The S75 epoch is only one transit, and it sits next to a gap.
- Insolation is 8.7 S⊕ (Teq about 437 K from the DV report). That is far inside the habitable zone; the "hz" in the lane name is a search label.
