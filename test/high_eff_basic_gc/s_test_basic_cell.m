% authors: bohan
% A simple script for testing the basic grating cell class
% Can be used to debug the grating cell class


clear; close all;

% Dependencies
% addpath( genpath( ['..' filesep '..' ] ));

% initial settings
disc            = 10;
units           = 'nm';
lambda          = 1310; %1500;
index_clad      = 1.45;
y_domain        = 2500;
fill_top        = 0.375;
fill_bot        = 0.55;
period          = 500;
offset_ratio    = 20.7/period;
layer_thick     = 100;

% make grating cell
GC = f_makeGratingCell_basic( disc, units, lambda, ...
                              index_clad, y_domain, period, ...
                              fill_top, fill_bot, ...
                              offset_ratio, layer_thick );
                                                          
% DEBUG plot the index
GC.plotIndex();

% run simulation
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
guessk      = 0.010217 + 0.000436i;
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [1, 100, 20, 2];

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

% Plot the accepted mode
figure;
imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (abs) for accepted mode, ka/2pi real = %f', real( GC.k*period/(2*pi) ) ) );

% display calculated k
fprintf('\nComplex k = %f + %fi\n', real(GC.k), imag(GC.k) );

% display radiated power
fprintf('\nRad power up = %e\n', GC.P_rad_up);
fprintf('Rad power down = %e\n', GC.P_rad_down);
fprintf('Up/down power directivity = %f\n', GC.directivity);

% display angle of radiation
fprintf('\nAngle of maximum radiation = %f deg\n', GC.max_angle_up);

% plot full Ez with grating geometry overlaid
GC.plot_E_field_gui();
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



























