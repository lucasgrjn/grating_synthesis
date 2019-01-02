#!/bin/bash -l

# Runs synthesis pipeline
# Inputs are lambda, optimal angle, and box thickness
#
# Example
#	qsub run_synth_bash.sh 1200 10 150
#		The above runs the synthesis for wavelength of 1200nm, angle 10 deg, and BOX thickness of 150nm

# set job name
#$ -N grating_synth

# set # cores to use
#$ -pe omp 28

# specify hard time limit for job
#$ -l h_rt=64:00:00

# merge output and error files in one
#$ -j y

# send me email if job finishes or aborted
#$ -m ea


# load matlab
module load matlab/2018a

# run script
matlab -nodisplay -r "addpath('/project/siphot/bz/code/utility'); cd(['..' filesep '45RFSOI']); f_run_synth_grating_v3( $1, $2, $3); exit"
