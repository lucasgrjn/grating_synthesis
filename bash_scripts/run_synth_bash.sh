#!/bin/bash -l

# set # cores to use
#$ -pe omp 28

# specify hard time limit for job
#$ -l h_rt=12:00:00

# merch output and error files in one
#$ -j y

# send me email if job finishes or aborted
#$ -m ea

# load matlab
module load matlab/2017b

# run script
matlab -nodisplay -r "f_run_param_sweep($NSLOTS); exit"
