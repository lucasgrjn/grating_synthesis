% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver

clear; close all;

% path to main code
addpath([ '..' filesep 'main']);

% set variables
N           = ones( 4, 4 );
disc        = 10;
k0          = 5;
num_modes   = 5;
guess_k     = 5;
BC          = 0;
PML_options = [ 0 0 0 0 ];

% run modesolver
[Phi_1D, k] = complexk_mode_solver_2D_PML_v2( N, disc, k0, num_modes, guess_k, BC, PML_options )

