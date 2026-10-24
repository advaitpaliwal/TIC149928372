#!/bin/bash
# usage: trc_launch_spark.sh [molusc] ; 10 runs (seeds 1000-1009) on the Spark, env as v4/vet_264151965/glue/trc_launch_spark.sh, nice 19.
cd ~/planet-hunt/trc_149928372; D=trc/out; [ "$1" = molusc ] && D=trc_molusc/out; mkdir -p $D
for k in 0 1 2 3 4 5 6 7 8 9; do
  NUMPY_MADVISE_HUGEPAGE=0 NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 nohup nice -n 19 ~/planet-hunt/trc_venv/bin/python -I glue/trc_run.py $k $1 > $D/log_$k.txt 2>&1 &
done
sleep 1; ps -eo args | grep -c "glue/trc_run.py"; date
