% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver

clear; close all;

% path to main code
addpath([ '..' filesep 'main']);

% set variables
N           = ones( 5, 5 );
disc        = 10;
k0          = 5;
num_modes   = 5;
guess_k     = 5;
BC          = 0;

% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
PML_options = [ 0 200 500 5 ];

% run modesolver
[Phi_1D, k] = complexk_mode_solver_2D_PML_v2( N, disc, k0, num_modes, guess_k, BC, PML_options )




% authors: bohan

% Script for testing and debugging the new two level grating cell class

clear; close all;

addpath(['..' filesep 'main']);

% initial settings
disc        = 5;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
domain      = [ 2000, 900 ];
numcells    = 10;

% Init a new object
GC = c_twoLevelGratingCell(  'discretization', disc, ...
                            'units', units, ...
                            'lambda', lambda, ...
                            'domain_size', domain, ...
                            'background_index', index_clad, ...
                            'numcells', numcells )

% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
fill            = 0.7;
ratio           = 1.0;
offset          = 0.0;
period          = domain(2);
wg_index        = [ 3.4, 3.4 ];
wg_thick        = [ 100, 100 ];
wg_min_y        = [ domain(1)/2, domain(1)/2-wg_thick(1) ];
wgs_duty_cycles = [ fill*ratio, fill ];
% wgs_offsets     = [ 0, offset*period ];
wgs_offsets     = [ 0, 200 ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, wg_index, ...
                                        wgs_duty_cycles, wgs_offsets );
                                 
                                 
% DEBUG plot the index
GC.plotIndex();

% run simulation
num_modes   = 1;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 500, 2 ];

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options );

% % DEBUG plot physical fields and all fields
% k_all       = GC.debug.k_all;
% ka_2pi_all  = k_all*domain(2)/(2*pi);
% 
% for ii = 1:length( k_all )
%     % Plotting physical fields
%     % plot field, abs
%     figure;
%     imagesc( GC.x_coords, GC.y_coords, abs( GC.debug.phi_all(:,:,ii) ) );
%     colorbar;
%     set( gca, 'YDir', 'normal' );
%     title( sprintf( 'Field (abs) for mode %i, k*a/2pi real = %f', ii, real( ka_2pi_all(ii) ) ) );
% end
% 
% % plot real and imag k
% k_labels = {};
% for ii = 1:length( k_all )
%     k_labels{end+1} = [ ' ', num2str(ii) ];
% end
% figure;
% plot( real( ka_2pi_all ), imag( ka_2pi_all ), 'o' ); 
% text( real( ka_2pi_all ), imag( ka_2pi_all ), k_labels );
% xlabel('real k*a/2pi'); ylabel('imag k*a/2pi');
% title('real vs imag k');
% makeFigureNice();

% Plot the accepted mode
figure;
imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (abs) for accepted mode, ka/2pi real = %f', real( GC.k*domain(2)/(2*pi) ) ) );

% display calculated k
fprintf('\nComplex k = %f + %fi\n', real(GC.k), imag(GC.k) );

% display radiated power
fprintf('\nRad power up = %e\n', GC.P_rad_up);
fprintf('Rad power down = %e\n', GC.P_rad_down);
fprintf('Up/down power directivity = %f\n', GC.directivity);

% display angle of radiation
fprintf('\nAngle of maximum radiation = %f deg\n', GC.max_angle_up);

% plot full Ez with grating geometry overlaid
GC.plotEz_w_edges();
axis equal;

% % plot a slice of Eup vs. Sup
% y_up    = 279;
% E_slice = GC.E_z( y_up, : );
% figure;
% plot( 1:length(E_slice), abs(E_slice)./max(abs(E_slice(:))) ); hold on;
% % plot( 1:length(E_slice), GC.debug.Sy_up./max(abs(GC.debug.Sy_up(:))) );
% plot( 1:length(E_slice), GC.debug.Sy_down./max(abs(GC.debug.Sy_down(:))) );
% legend('E field (abs)', 'S_y');
% title('E field and S_y, normalized to 1');
% makeFigureNice();
% 
% % plot input S and E
% E_in = GC.E_z( :, 1 );
% figure;
% plot( 1:length(E_in), abs(E_in)./max(abs(E_in(:))) ); hold on;
% plot( 1:length(E_in), GC.debug.Sx_in./max(abs(GC.debug.Sx_in(:))) );
% % plot( 1:length(E_slice), GC.debug.Sy_down./max(abs(GC.debug.Sy_down(:))) );
% legend('E field (abs)', 'S_x');
% title('E field and S_x at input, normalized to 1');
% makeFigureNice();
% 



























