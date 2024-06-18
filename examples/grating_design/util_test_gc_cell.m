% This script shows how you can test your grating cell

clear; close all;

% grating synthesis codes
addpath(genpath('../../..'));

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

% guess what the neff and period should be
neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 15; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;

dutycycle = 0.9;

GC = f_makeGratingCell( ...
    dxy, period, dutycycle );
                            
% plot index
GC.plotIndex();

% set grating solver settings
num_modes   = 15;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2]; % [ <on or off>, thickness, strength, polynomial order ]
guessk      = k0 * neff_slab;
OPTS        = struct(); % optional args

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );

% plot e field
GC.plot_E_field_gui();