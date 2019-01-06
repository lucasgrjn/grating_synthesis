#!/bin/bash -l

# Runs synthesis pipeline for a basic grating coupler
# Inputs are lambda, optimal angle, layer thickness
#
# Example
#	qsub run_synth_bash_basic_gc_topbot.sh 1200 -10 100 
#		The above runs the synthesis for wavelength of 1200nm, angle -10 and layer thickness of 100nm

# set job name
#$ -N grating_synth_topbot

# set # cores to use
#$ -pe omp 28

# specify hard time limit for job
#$ -l h_rt=23:00:00

# merge output and error files in one
#$ -j y

# send me email if job finishes or aborted
#$ -m ea


# load matlab
module load matlab/2018b

# run script
matlab -nodisplay -r "addpath('/project/siphot/bz/code/utility'); addpath( genpath ( ['/project/siphot/bz/git/grating_synthesis'] ) ); f_run_synth_grating_basic_gc_topbot( $1, $2, $3 ); exit"
