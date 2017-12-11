% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% debugging by comparing with analytical bragg stack result

clear; close all;

% Script for testing and debugging the new two level grating cell class

clear; close all;

% import code
addpath(['..' filesep 'main']);         % main
% addpath(['..' filesep '45RFSOI']);      % 45rfsoi

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1000; %1500;
index_clad  = 1.0;

% make index
% let's have a quarter wave stack at 1000nm
% indices
n1      = 1.0;
n2      = 1.25;
d1      = lambda/(n1*4);
d2      = lambda/(n2*4);
period  = d1+d2;
domain  = [ 100, period ];

% draw indices
% x = dir of propagation, y = transverse
x_coords    = 0:disc:period-disc;
y_coords    = 0:disc:domain(1)-disc;
N           = zeros( length(y_coords), length(x_coords) );
N( :, x_coords >= 0 & x_coords < d1 )   = n1;
N( :, x_coords >= d1 )                  = n2;

% DEBUG plot N
figure;
imagesc( x_coords, y_coords, N );
xlabel('x'); ylabel('y');
set( gca, 'ydir', 'normal' );
colorbar;
title('DEBUG N');
                                 

% run simulation
num_modes   = 20;
BC          = 1;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 0, 200, 500, 2 ];


% solve bandstructure
lambda_min  = 200;
lambda_max  = 5000;
k0_all      = linspace( 2*pi/lambda_min, 2*pi/lambda_max, 10 );

% run solver
guessk = pi/(2*period);
for ii = 1:length(k0_all)
    
    fprintf('\nloop %i of %i\n\n', ii, length(k0_all));
    
    % run old
    fprintf('running old solver\n');
    tic;
    [Phi_1D_old, k_old, Phi_all_old, k_all_old, A_old, B_old] = complexk_mode_solver_2D_PML_old( N, ...
                                                                                   disc, ...
                                                                                   k0_all(ii), ...
                                                                                   num_modes, ...
                                                                                   guessk, ...
                                                                                   BC, ...
                                                                                   pml_options );
    toc;

    % run new
    fprintf('running new solver\n');
    tic;
    [Phi_all, k_all, A, B] = complexk_mode_solver_2D_PML( N, ...
                                                           disc, ...
                                                           k0_all(ii), ...
                                                           num_modes, ...
                                                           guessk, ...
                                                           BC, ...
                                                           pml_options );
    toc;
    
end


                                                                           
% % reshape and sort the Phis
% % when they come out raw from the modesolver, Phi_all's columns are the
% % eigenvectors
% % The eigenvectors are wrapped by column, then row
% ny = domain(1)/disc;
% nx = domain(2)/disc;
% Phi_all_half    = Phi_all( 1:end/2, : );                                        % first remove redundant bottom half
% Phi_all_reshape = reshape( Phi_all_half, ny, nx, size(Phi_all_half, 2) );       % hopefully this is dimensions y vs. x vs. mode#
% Phi_firstmode   = Phi_all_reshape( :, :, 1 );
% Phi_secondmode  = Phi_all_reshape( :, :, 2 );

% do the same but with the old phi for comparison
Phi_all_half_old    = Phi_all_old( 1:end/2, : );
Phi_all_old_reshape = reshape( Phi_all_half_old, nx, ny, size(Phi_all_half_old, 2) );       % hopefully this is dimensions x vs. y vs. mode#
Phi_firstmod_old    = Phi_all_old_reshape( :, :, 1 ).';                                     % gotta transpose

% x and y coords
x_coords = 0:disc:domain(2)-disc;
y_coords = 0:disc:domain(1)-disc;

% DEBUG plot firstmode
figure;
imagesc( x_coords, y_coords, real( Phi_firstmode ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (real) for mode 1, new ver' ) );
% DEBUG plot firstmode
figure;
imagesc( x_coords, y_coords, real( Phi_firstmod_old ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (real) for mode 1, old ver' ) );

% DEBUG plot secondmode
figure;
imagesc( x_coords, y_coords, real( Phi_secondmode ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (real) for mode 2, new ver' ) );

% debug display imag(k)
imag(k_all)

% % plot the sparse distributions
% figure;
% spy(A);
% title('A matrix, new');
% 
% figure;
% spy(A_old);
% title('A matrix, old');
% 
% figure;
% spy(B);
% title('B matrix, new');
% 
% figure;
% spy(B_old);
% title('B matrix, old');

% % re-scale k
% k = k * nm * obj.units.scale;
% 
% % run simulation
% GC = GC.runSimulation( num_modes, BC, pml_options );

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

% % Plot the accepted mode
% figure;
% imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi ) );
% colorbar;
% set( gca, 'YDir', 'normal' );
% title( sprintf( 'Field (abs) for accepted mode, ka/2pi real = %f', real( GC.k*domain(2)/(2*pi) ) ) );
% 
% % display calculated k
% fprintf('\nComplex k = %f + %fi\n', real(GC.k), imag(GC.k) );
% 
% % display radiated power
% fprintf('\nRad power up = %e\n', GC.P_rad_up);
% fprintf('Rad power down = %e\n', GC.P_rad_down);
% fprintf('Up/down power directivity = %f\n', GC.directivity);
% 
% % display angle of radiation
% fprintf('\nAngle of maximum radiation = %f deg\n', GC.max_angle_up);
% 
% % plot full Ez with grating geometry overlaid
% GC.plotEz_w_edges();
% axis equal;

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



























