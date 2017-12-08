#!/bin/bash -l

# OLD run a bunch of sweeps
# qsub -v fill1=0.02 -v fill2=0.04 run_synth_bash_fill.sh
# qsub -v fill1=0.06 -v fill2=0.08 run_synth_bash_fill.sh
# qsub -v fill1=0.10 -v fill2=0.12 run_synth_bash_fill.sh
# qsub -v fill1=0.14 -v fill2=0.16 run_synth_bash_fill.sh
# qsub -v fill1=0.18 -v fill2=0.20 run_synth_bash_fill.sh
# qsub -v fill1=0.22 -v fill2=0.24 run_synth_bash_fill.sh
# qsub -v fill1=0.26 -v fill2=0.28 run_synth_bash_fill.sh
# qsub -v fill1=0.30 -v fill2=0.32 run_synth_bash_fill.sh
# qsub -v fill1=0.34 -v fill2=0.36 run_synth_bash_fill.sh
# qsub -v fill1=0.38 -v fill2=0.40 run_synth_bash_fill.sh
# qsub -v fill1=0.42 -v fill2=0.44 run_synth_bash_fill.sh
# qsub -v fill1=0.46 -v fill2=0.48 run_synth_bash_fill.sh
# qsub -v fill1=0.50 -v fill2=0.52 run_synth_bash_fill.sh
# qsub -v fill1=0.54 -v fill2=0.56 run_synth_bash_fill.sh
# qsub -v fill1=0.58 -v fill2=0.60 run_synth_bash_fill.sh
# qsub -v fill1=0.62 -v fill2=0.64 run_synth_bash_fill.sh
# qsub -v fill1=0.66 -v fill2=0.68 run_synth_bash_fill.sh
# qsub -v fill1=0.70 -v fill2=0.72 run_synth_bash_fill.sh
# qsub -v fill1=0.74 -v fill2=0.76 run_synth_bash_fill.sh
# qsub -v fill1=0.78 -v fill2=0.80 run_synth_bash_fill.sh
# qsub -v fill1=0.82 -v fill2=0.84 run_synth_bash_fill.sh
# qsub -v fill1=0.86 -v fill2=0.88 run_synth_bash_fill.sh
# qsub -v fill1=0.90 -v fill2=0.92 run_synth_bash_fill.sh
# qsub -v fill1=0.94 -v fill2=0.96 run_synth_bash_fill.sh
# qsub -v fill1=0.98 -v fill2=0.99 run_synth_bash_fill.sh
# qsub -v fill1=0.98 -v fill2=0.99 run_synth_bash_fill.sh
# qsub -v fill1=0.98 -v fill2=0.99 run_synth_bash_fill.sh
# qsub -v fill1=0.98 -v fill2=0.99 run_synth_bash_fill.sh

# NEW run a bunch of sweeps, one fill at a time
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

# qsub -v fill=0.30 run_synth_bash_fill.sh
# qsub -v fill=0.31 run_synth_bash_fill.sh
# qsub -v fill=0.32 run_synth_bash_fill.sh
# qsub -v fill=0.33 run_synth_bash_fill.sh
# qsub -v fill=0.34 run_synth_bash_fill.sh
# qsub -v fill=0.35 run_synth_bash_fill.sh
qsub -v fill=0.36 run_synth_bash_fill.sh
# qsub -v fill=0.37 run_synth_bash_fill.sh
# qsub -v fill=0.38 run_synth_bash_fill.sh
# qsub -v fill=0.39 run_synth_bash_fill.sh

qsub -v fill=0.40 run_synth_bash_fill.sh
# qsub -v fill=0.41 run_synth_bash_fill.sh
qsub -v fill=0.42 run_synth_bash_fill.sh
# qsub -v fill=0.43 run_synth_bash_fill.sh
# qsub -v fill=0.44 run_synth_bash_fill.sh
# qsub -v fill=0.45 run_synth_bash_fill.sh
qsub -v fill=0.46 run_synth_bash_fill.sh
# qsub -v fill=0.47 run_synth_bash_fill.sh
qsub -v fill=0.48 run_synth_bash_fill.sh
# qsub -v fill=0.49 run_synth_bash_fill.sh

qsub -v fill=0.50 run_synth_bash_fill.sh
# qsub -v fill=0.51 run_synth_bash_fill.sh
qsub -v fill=0.52 run_synth_bash_fill.sh
# qsub -v fill=0.53 run_synth_bash_fill.sh
qsub -v fill=0.54 run_synth_bash_fill.sh
# qsub -v fill=0.55 run_synth_bash_fill.sh
qsub -v fill=0.56 run_synth_bash_fill.sh
# qsub -v fill=0.57 run_synth_bash_fill.sh
qsub -v fill=0.58 run_synth_bash_fill.sh
# qsub -v fill=0.59 run_synth_bash_fill.sh

qsub -v fill=0.60 run_synth_bash_fill.sh
# qsub -v fill=0.61 run_synth_bash_fill.sh
qsub -v fill=0.62 run_synth_bash_fill.sh
# qsub -v fill=0.63 run_synth_bash_fill.sh
qsub -v fill=0.64 run_synth_bash_fill.sh
# qsub -v fill=0.65 run_synth_bash_fill.sh
# qsub -v fill=0.66 run_synth_bash_fill.sh
# qsub -v fill=0.67 run_synth_bash_fill.sh
qsub -v fill=0.68 run_synth_bash_fill.sh
# qsub -v fill=0.69 run_synth_bash_fill.sh

# qsub -v fill=0.70 run_synth_bash_fill.sh
# qsub -v fill=0.71 run_synth_bash_fill.sh
# qsub -v fill=0.72 run_synth_bash_fill.sh
# qsub -v fill=0.73 run_synth_bash_fill.sh
# qsub -v fill=0.74 run_synth_bash_fill.sh
# qsub -v fill=0.75 run_synth_bash_fill.sh
# qsub -v fill=0.76 run_synth_bash_fill.sh
# qsub -v fill=0.77 run_synth_bash_fill.sh
qsub -v fill=0.78 run_synth_bash_fill.sh
# qsub -v fill=0.79 run_synth_bash_fill.sh

# qsub -v fill=0.80 run_synth_bash_fill.sh
# qsub -v fill=0.81 run_synth_bash_fill.sh
# qsub -v fill=0.82 run_synth_bash_fill.sh
# qsub -v fill=0.83 run_synth_bash_fill.sh
# qsub -v fill=0.84 run_synth_bash_fill.sh
# qsub -v fill=0.85 run_synth_bash_fill.sh
# qsub -v fill=0.86 run_synth_bash_fill.sh
# qsub -v fill=0.87 run_synth_bash_fill.sh
# qsub -v fill=0.88 run_synth_bash_fill.sh
# qsub -v fill=0.89 run_synth_bash_fill.sh

# qsub -v fill=0.90 run_synth_bash_fill.sh
# qsub -v fill=0.91 run_synth_bash_fill.sh
# qsub -v fill=0.92 run_synth_bash_fill.sh
# qsub -v fill=0.93 run_synth_bash_fill.sh
# qsub -v fill=0.94 run_synth_bash_fill.sh
# qsub -v fill=0.95 run_synth_bash_fill.sh
# qsub -v fill=0.96 run_synth_bash_fill.sh
# qsub -v fill=0.97 run_synth_bash_fill.sh
qsub -v fill=0.98 run_synth_bash_fill.sh
# qsub -v fill=0.99 run_synth_bash_fill.sh

# move all logfiles
# mv ./run_synth_bash_fill.sh.o* ./logfiles/run_synth_bash_fill


