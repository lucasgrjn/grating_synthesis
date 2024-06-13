% testing the simplest possible grating with TM

clear; close all;

% import code
addpath( genpath( '..' ) );

% slab mode solver
addpath('G:\My Drive\research\popovic group\code\slab_modesolver');

% initial settings (nm)
disc        = 10;
lambda      = 1550;
k0          = 2*pi/lambda;
index_clad  = 1.45;
index_core  = 3.47;
y_size      = 4000;

% grating dimensions
si_thick    = 200;
gap_len     = 250;

% first get slab mode k
[ neff, beta, kx, ~, ~ ] = solve_symm_slab( si_thick, index_clad, index_core, ...
                                            lambda, 'TM', false );

% solve for period for desired angle
theta   = 15 * pi/180;
period  = (2*pi) ./ ( beta(1) - k0*sin(theta) );
                                        
% draw index
z_coords    = 0 : disc : period-disc;
x_coords    = 0 : disc : y_size-disc;
N           = index_clad * ones( length(x_coords), length(z_coords) );
x_high      = x_coords > x_coords(end/2) - si_thick/2 - disc/2;             % waveguide
x_low       = x_coords < x_coords(end/2) + si_thick/2 - disc/2;             % waveguide
% x_right     = x_coords > gap_len - disc/2;                                  % start with the gap on the left
z_left     = z_coords < period - gap_len - disc/2;                                  
N( x_high & x_low, z_left ) = index_core;

% plot index
figure;
imagesc( z_coords, x_coords, N );
xlabel('z (nm)'); ylabel('x (nm)');
set( gca, 'ydir', 'normal' );
colorbar;
axis image;
title('index distribution');
          

% modesolver settings
guessk      = beta(1);
guessk      = 0.0069
num_modes   = 10;
BC          = 1;        % 0 for PEC, 1 for PMC
pol         = 'TM';     % 'TE' or 'TM'
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 10, 2 ];

% run solver
tic;
[Phi_solved, k_solved, ~, ~] = f_bloch_complexk_mode_solver_2D_PML( N, ...
                                                       disc, ...
                                                       k0, ...
                                                       num_modes, ...
                                                       guessk, ...
                                                       BC, ...
                                                       pol, ...
                                                       pml_options );
toc;


f_plot_all_modes_gui(  Phi_solved, z_coords, x_coords, k_solved );

















