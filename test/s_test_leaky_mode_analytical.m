% script for testing leaky mode solution function
%
% authors: bohan zhang

clear; close all;

% add path to solver
addpath( ['..' filesep 'auxiliary_functions'] );

% define inputs
lambda0 = 1000;
a       = 0.39*lambda0;
b       = 1.96*lambda0;
n0      = 1.45;
n1      = 1.39;


% run solver
f_solve_leaky_modes_analytical( a, b, n0, n1, lambda0 )