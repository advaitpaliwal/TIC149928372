#!/bin/zsh
# 4 MOLUSC runs (serial, 1 thread) (seeds 11-14), 50,000 companions each, RUWE + Gaia contrast only, age default 5 Gyr; same command form as TIC 240488148.
cd ~/Scratch/planet-hunt/hz-bright/vet_149928372
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
P=~/Scratch/planet-hunt/tic240488148-companions/.venv-molusc/bin/python
O=$PWD/molusc_out
run(){ SEED=$1 taskpolicy -a $P glue/molusc_run.py --ruwe --gaia -a $O/m50k_s$1 10h48m16.058s 65d06m00.38s 50000 0.62 > $O/m50k_s$1.log 2>&1 }
run 11; run 12; run 13; run 14
echo DONE > $O/DONE
