#!/bin/zsh
# Fetch the 7 light-curve products listed in pool.csv (SPOC 2-min S14,S21; TESS-SPOC FFI S41,47,48,74,75) from MAST.
cd ~/Scratch/planet-hunt/hz-bright/vet_149928372
B=https://mast.stsci.edu/api/v0.1/Download/file?uri=
for pair in 14:mast:TESS/product/tess2019198215352-s0014-0000000149928372-0150-s_lc.fits 21:mast:TESS/product/tess2020020091053-s0021-0000000149928372-0167-s_lc.fits; do
 s=${pair%%:*}; u=$B${pair#*:}; curl -sfL -m 120 --retry 3 -o lcfits/spoc2m_s$s.fits "$u" && echo "$u" > lcfits/spoc2m_s$s.fits.url; done
for s in 41 47 48 74 75; do S=$(printf %04d $s)
 u=https://archive.stsci.edu/hlsps/tess-spoc/s$S/target/0000/0001/4992/8372/hlsp_tess-spoc_tess_phot_0000000149928372-s${S}_tess_v1_lc.fits
 curl -sfL -m 120 --retry 3 -o lcfits/spoc_s$s.fits $u && echo "$u" > lcfits/spoc_s$s.fits.url; done
ls -la lcfits
