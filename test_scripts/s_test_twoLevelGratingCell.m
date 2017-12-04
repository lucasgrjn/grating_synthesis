% authors: bohan

% Script for testing and debugging the new two level grating cell class

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
% domain      = [ 1600, 470 ];
domain      = [ 1600, 700 ];

% % quick back of envelope calculation for period
% theta   = 15 * pi/180;
% neff    = 2.25;    % very approximate
% kg      = 2*pi*neff/lambda;
% kclad   = 2*pi*index_clad/lambda;
% period_approx = (2*pi)/( kg - kclad*sin(theta) );
% 
% % round the period and set the domain
% period_approx_round = round( period_approx * 1e-1 )/1e-1;
% period              = period_approx_round;
% domain              = [ 1600, period ];

% Init a new object
Q = c_twoLevelGratingCell(  'discretization', disc, ...
                            'units', units, ...
                            'lambda', lambda, ...
                            'domain_size', domain, ...
                            'background_index', index_clad )


% draw two levels using two level builder function
wg_thick        = [ 100, 100 ];
wg_min_y        = [ domain(1)/2-wg_thick(1), domain(1)/2 ];
wg_indx         = [ 3.4, 3.4 ];
wgs_duty_cycles = [ 0.5, 0.5 ];
wgs_offsets     = [ 0, 600 ];
Q               = Q.twoLevelBuilder( wg_min_y, wg_thick, wg_indx, ...
                                     wgs_duty_cycles, wgs_offsets );

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

% % % DEBUG display the unguided power
% % Q.debug.unguided_power
% Q.debug.P_rad_up_Hx_only
% Q.debug.P_rad_up_Hx_Hy 
% Q.debug.P_rad_down_Hx_only
% Q.debug.P_rad_down_Hx_Hy

% display calculated k
fprintf('\nComplex k = %f + %fi\n', real(Q.k), imag(Q.k) );

% display radiated power
fprintf('\nRad power up = %e\n', Q.P_rad_up);
fprintf('Rad power down = %e\n', Q.P_rad_down);
fprintf('Up/down power directivity = %f\n', Q.directivity);

% display angle of radiation
fprintf('\nAngle of maximum radiation = %f deg\n', Q.max_angle_up);







