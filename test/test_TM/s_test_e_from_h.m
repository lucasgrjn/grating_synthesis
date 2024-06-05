% Testing solving E field from H field
% using kenaish's phc example

clear; close all;

% import code
addpath( genpath( '..' ) );

% initial settings (nm)
disc        = 10;
lambda      = 1000; %1500;
index_clad  = 1.45;
index_core  = 3.45;
y_size      = 8000;

% PhC dimensions
period      = 450;
si_width    = 600;
gap_len     = 220;

% draw index
x_coords    = 0:disc:period-disc;
y_coords    = 0:disc:y_size-disc;
N           = index_clad * ones( length(y_coords), length(x_coords) );
y_high      = y_coords > y_coords(end/2) - si_width/2 - disc/2;             % waveguide
y_low       = y_coords < y_coords(end/2) + si_width/2 - disc/2;             % waveguide
x_right     = x_coords > gap_len - disc/2;                                  % start with the gap on the left
N( y_high & y_low, x_right ) = index_core;

% plot index
figure;
imagesc( x_coords, y_coords, N );
xlabel('x (nm)'); ylabel('y (nm)');
set( gca, 'ydir', 'normal' );
colorbar;
% axis image;
title('index distribution');
          

% modesolver settings
num_modes   = 1;
BC          = 0;        % 0 for PEC, 1 for PMC
pol         = 'TM';     % 'TE' or 'TM'
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 20, 2 ];
% k0          = 0.0021;    
k0          = 0.35 * pi/period;
% guessk      = 0.005;
guessk      = 1 * pi/period;

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
    

% plot mode in gui
f_plot_all_modes_gui( Phi_solved, x_coords, y_coords, k_solved );
    
% compute Hz, by adding in the phase
Hz   = Phi_solved .* repmat( exp( 1i * k_solved * x_coords ), length(y_coords), 1 );

% DEBUG plot Hz
figure;
imagesc( x_coords, y_coords, real(Hz) );
xlabel('x'); ylabel('y');
colorbar;
colormap('redbluehilight');
set(gca, 'ydir', 'normal');
title('H_z (real)');

% compute Ex, Ey
[ Ex, Ey ] = f_get_E_from_H_field( Hz, N );
    
% half grid x and y
x_coords_halfgrid = (x_coords(2:end) + x_coords(1:end-1))/2;
y_coords_halfgrid = (y_coords(2:end) + y_coords(1:end-1))/2;

% plot Ex
figure;
imagesc( x_coords_halfgrid, y_coords_halfgrid, real(Ex) );
xlabel('x'); ylabel('y');
colorbar;
colormap('redbluehilight');
set(gca, 'ydir', 'normal');
title('E_x (real)');

% plot Ey
figure;
imagesc( x_coords_halfgrid, y_coords_halfgrid, real(Ey) );
xlabel('x'); ylabel('y');
colorbar;
colormap('redbluehilight');
set(gca, 'ydir', 'normal');
title('E_y (real)');
























