#!/bin/bash -l

# set # cores to use
#$ -pe omp 16

# specify hard time limit for job
#$ -l h_rt=12:00:00

# merch output and error files in one
#$ -j y

# send me email if job finishes or aborted
#$ -m ea

# pick fills
#$ -v fill=0.01

# load matlab
module load matlab/2017b

# run script
matlab -nodisplay -r "addpath(['..' filesep 'main']); addpath(['..' filesep '45RFSOI']); f_run_param_sweep_45RFSOI_debug($NSLOTS, $fill); exit"
