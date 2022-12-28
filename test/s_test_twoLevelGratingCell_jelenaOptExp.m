% authors: bohan

% Script for testing and debugging the new two level grating cell class

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
domain      = [ 1600, 470 ];
% domain      = [ 1600, 700 ];

% quick back of envelope calculation for period
theta   = 80 * pi/180;
neff    = 2.8;    % very approximate
kg      = 2*pi*neff/lambda;
kclad   = 2*pi*index_clad/lambda;
period_approx = (2*pi)/( kg - kclad*sin(theta) );


% Init a new object
Q = c_twoLevelGratingCell(  'discretization', disc, ...
                            'units', units, ...
                            'lambda', lambda, ...
                            'domain_size', domain, ...
                            'background_index', index_clad )

% Add a layer
height_y    = 260;
min_y       = (domain(1)-height_y)/2;
index       = 3.4;
Q           = Q.addLayer( min_y, height_y, index );
Q.wg_min_y  = min_y;
Q.wg_max_y  = min_y + height_y;

% add first rectangle
width_x     = 100;
min_x       = 0;
min_y       = min_y+20;
height_y    = 240;
index       = index_clad;
Q           = Q.addRect( min_x, min_y, width_x, height_y, index );

% add second rectangle
width_x     = 150;
min_x       = 130;
min_y       = min_y + (240-60);
height_y    = 60;
index       = index_clad;
Q           = Q.addRect( min_x, min_y, width_x, height_y, index );

% DEBUG plot the index
Q.plotIndex();

% run simulation
num_modes   = 20;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 500, 2 ];

% run simulation
Q = Q.runSimulation( num_modes, BC, pml_options );

% DEBUG plot physical fields and all fields
k_all       = Q.debug.k_all;
ka_2pi_all  = k_all*domain(2)/(2*pi);

for ii = 1:length( k_all )
    % Plotting physical fields
    % plot field, abs
    figure;
    imagesc( Q.x_coords, Q.y_coords, abs( Q.debug.phi_all(:,:,ii) ) );
    colorbar;
    set( gca, 'YDir', 'normal' );
    title( sprintf( 'Field (abs) for mode %i, k*a/2pi real = %f', ii, real( ka_2pi_all(ii) ) ) );
end

% plot real and imag k
k_labels = {};
for ii = 1:length( k_all )
    k_labels{end+1} = [ ' ', num2str(ii) ];
end
figure;
plot( real( ka_2pi_all ), imag( ka_2pi_all ), 'o' ); 
text( real( ka_2pi_all ), imag( ka_2pi_all ), k_labels );
xlabel('real k*a/2pi'); ylabel('imag k*a/2pi');
title('real vs imag k');
makeFigureNice();

% Plot the accepted mode
figure;
imagesc( Q.x_coords, Q.y_coords, abs( Q.Phi ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (abs) for accepted mode, ka/2pi real = %f', real( Q.k*domain(2)/(2*pi) ) ) );

% % DEBUG display the unguided power
% Q.debug.unguided_power
Q.debug.P_rad_up_Hx_only
Q.debug.P_rad_up_Hx_Hy 
Q.debug.P_rad_down_Hx_only
Q.debug.P_rad_down_Hx_Hy






