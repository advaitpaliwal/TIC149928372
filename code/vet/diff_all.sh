#!/bin/zsh
# usage: diff_all.sh SECTOR... ; 3 at a time, single-threaded
cd ~/Scratch/planet-hunt/hz-bright/vet_149928372
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
for s in "$@"; do echo $s; done | xargs -P 3 -I{} sh -c 'cd glue/.. && taskpolicy -a ~/Scratch/planet-hunt/hunt4/.venv/bin/python -I -c "import sys; sys.path.insert(0,\"glue\"); sys.argv=[\"diff.py\",\"{}\"]; import diff; diff.main()" >> pix/diff.log 2>&1'
