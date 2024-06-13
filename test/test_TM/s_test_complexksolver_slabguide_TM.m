% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% for the TM mode

clear; close all;

% dependencies
% bloch solver
addpath( genpath( '..' ) );
% slab solver
addpath( 'C:\Users\beezy\Google Drive\research\popovic group\code\slab_modesolver' );
addpath( 'D:\Google Drive\research\popovic group\code\slab_modesolver' );
addpath( 'G:\My Drive\research\popovic group\code\slab_modesolver' );

% initial settings
disc        = 2.5;
units       = 'nm';
lambda      = 1000; %1500;
index_clad  = 1.0;
k0          = 2*pi/lambda;
pol         = 'TM';

% make index
% would like to make a single mode wg
n1      = 1.5;
n2      = 2.0;                  % 1.25;
t_wg    = 500;
period  = disc*10;
domain  = [ 2000, period ];     % dimensions x, z

% draw indices
% % z = dir of propagation, x = transverse
z_coords    = 0:disc:domain(2)-disc;
x_coords    = 0:disc:domain(1)-disc;
N           = n1*ones( length(x_coords), length(z_coords) );
x_indx_wg   = ( x_coords > domain(1)/2 - t_wg/2 - disc/2 ) & (x_coords < domain(1)/2 + t_wg/2 - disc/2 );
N( x_indx_wg, : )   = n2;
% TESTING SOME INDEX AVERAGING
% N( 1:end-1, : ) = ( N( 1:end-1, : ) + N( 2:end, : ) )/2;
% N = imfilter( N, (1/9) * [ 1, 1, 1; 1, 1, 1; 1, 1, 1 ], 'replicate' );

% DEBUG plot N
figure;
imagesc( z_coords, x_coords, N );
xlabel('z'); ylabel('x');
set( gca, 'ydir', 'normal' );
colorbar;
title('DEBUG N');
        
% compute analytical solution (symmetric waveguide)
core_d      = t_wg;
n_clad      = n1;
n_core      = n2;
lambda0     = lambda;
[ neff, k_analytical, kx_temp, alpha_temp, field_analy ] = solve_symm_slab( core_d, n_clad, n_core, lambda0, pol, true );

% convert units to nm
% k_analytical = k_analytical;

% run simulation
% guessk      = k0*(n1+n2)/2;
guessk      = k_analytical(1);
num_modes   = 5;
BC          = 1;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 0, 200, 50, 2 ];

% run solver
fprintf('Running bloch solver\n');
tic;
[Phi_all, k_all, A, B] = f_bloch_complexk_mode_solver_2D_PML(  N, ...
                                                               disc, ...
                                                               k0, ...
                                                               num_modes, ...
                                                               guessk, ...
                                                               BC, ...
                                                               pol, ...
                                                               pml_options );
toc;


% plot mode
figure;
imagesc( z_coords, x_coords, real(Phi_all(:,:,1)) );
colorbar;
xlabel('x'); ylabel('y');
set( gca, 'ydir', 'normal' );
title('first mode field, real');
% plot mode
figure;
imagesc( z_coords, x_coords, abs(Phi_all(:,:,1)) );
colorbar;
xlabel('x'); ylabel('y');
set( gca, 'ydir', 'normal' );
title('first mode field, amp');


% plot all modes
f_plot_all_modes_gui(  Phi_all, z_coords, x_coords, k_all );


% Compare analytical with mode solved fields
figure;
% analytical field
plot( field_analy.x, field_analy.Hy(1,:) ); hold on;
% mode solved field
plot( x_coords - x_coords(end/2), abs(Phi_all( :, 1, 1 ))./max(abs(Phi_all( :, 1, 1 ))) );
% index distribution
plot( x_coords - x_coords(end/2), N(:,1) );
xlabel('x (nm)'); ylabel('H_y');
legend('analytical', 'mode solved', 'index');
title('Analytical vs. mode solver field profile, 1st mode');
makeFigureNice();

% % what if i add mode 1 and 2
% Phi_1_2 = Phi_all( :, :, 1 ) + Phi_all( :, :, 2);
% % plot mode
% figure;
% imagesc( x_coords, y_coords, real(Phi_1_2) );
% colorbar;
% xlabel('x'); ylabel('y');
% set( gca, 'ydir', 'normal' );
% title('first + 2nd mode field, real');
% % plot mode
% figure;
% imagesc( x_coords, y_coords, abs(Phi_1_2) );
% colorbar;
% xlabel('x'); ylabel('y');
% set( gca, 'ydir', 'normal' );
% title('first + 2nd mode field, amp');
% 
% % Compare analytical with mode solved fields
% figure;
% % analytical field
% plot( field_analy.x, field_analy.Hy(1,:) ); hold on;
% % mode solved field
% plot( y_coords - y_coords(end/2), abs(Phi_1_2( :, 1, 1 ))./max(abs(Phi_1_2( :, 1, 1 ))) );
% % index distribution
% plot( y_coords - y_coords(end/2), N(:,1) );
% xlabel('x (nm)'); ylabel('H_y');
% legend('analytical', 'mode solved', 'index');
% title('Analytical vs. mode solver field profile, 1st +2nd mode');
% makeFigureNice();



















