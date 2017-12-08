#!/bin/bash -l

# run a bunch of sweeps, one fill (top fill) at a time
# qsub -v fill=0.01 run_synth_bash_fill.sh
# qsub -v fill=0.02 run_synth_bash_fill.sh
# qsub -v fill=0.03 run_synth_bash_fill.sh
# qsub -v fill=0.04 run_synth_bash_fill.sh
# qsub -v fill=0.05 run_synth_bash_fill.sh
# qsub -v fill=0.06 run_synth_bash_fill.sh
# qsub -v fill=0.07 run_synth_bash_fill.sh
# qsub -v fill=0.08 run_synth_bash_fill.sh
# qsub -v fill=0.09 run_synth_bash_fill.sh

# qsub -v fill=0.10 run_synth_bash_fill.sh
# qsub -v fill=0.11 run_synth_bash_fill.sh
# qsub -v fill=0.12 run_synth_bash_fill.sh
# qsub -v fill=0.13 run_synth_bash_fill.sh
# qsub -v fill=0.14 run_synth_bash_fill.sh
# qsub -v fill=0.15 run_synth_bash_fill.sh
# qsub -v fill=0.16 run_synth_bash_fill.sh
# qsub -v fill=0.17 run_synth_bash_fill.sh
# qsub -v fill=0.18 run_synth_bash_fill.sh
# qsub -v fill=0.19 run_synth_bash_fill.sh

# qsub -v fill=0.20 run_synth_bash_fill.sh
# qsub -v fill=0.21 run_synth_bash_fill.sh
# qsub -v fill=0.22 run_synth_bash_fill.sh
# qsub -v fill=0.23 run_synth_bash_fill.sh
# qsub -v fill=0.24 run_synth_bash_fill.sh
# qsub -v fill=0.25 run_synth_bash_fill.sh
# qsub -v fill=0.26 run_synth_bash_fill.sh
# qsub -v fill=0.27 run_synth_bash_fill.sh
# qsub -v fill=0.28 run_synth_bash_fill.sh
# qsub -v fill=0.29 run_synth_bash_fill.sh

qsub -v fill=0.30 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.32 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.34 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.36 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.38 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.40 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.42 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.44 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.46 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.48 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.50 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.52 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.54 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.56 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.58 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.60 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.62 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.64 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.66 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.68 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.70 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.72 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.74 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.76 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.78 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.80 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.82 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.84 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.86 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.88 run_synth_bash_fill_45RFSOI.sh

qsub -v fill=0.90 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.92 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.94 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.96 run_synth_bash_fill_45RFSOI.sh
qsub -v fill=0.98 run_synth_bash_fill_45RFSOI.sh

# move all logfiles
# mv ./run_synth_bash_fill.sh.o* ./logfiles/run_synth_bash_fill