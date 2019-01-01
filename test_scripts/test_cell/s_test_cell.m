% authors: bohan
% A simple script for testing the basic grating cell class
% Can be used to debug the grating cell class

clear; close all;

% Dependencies
addpath( genpath( ['..' filesep '..' ] ));

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
domain      = [ 2000, 800 ];
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
fill            = 0.8;
ratio           = 1.0;
offset          = 0.0;
period          = domain(2);
wg_index        = [ 3.4, 3.4 ];
wg_thick        = [ 100, 100 ];
wg_min_y        = [ domain(1)/2, domain(1)/2-wg_thick(1) ];
wgs_duty_cycles = [ fill*ratio, fill ];
wgs_offsets     = [ 0, 200 ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, wg_index, ...
                                        wgs_duty_cycles, wgs_offsets );
                                 
                                 
% DEBUG plot the index
GC.plotIndex();

% run simulation
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
guessk      = 2*pi*wg_index(1)*fill/lambda;
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 5, 2 ];

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

% DEBUG plot physical fields and all fields
k_all       = GC.k_vs_mode;
ka_2pi_all  = k_all*domain(2)/(2*pi);

for ii = 1:length( k_all )
    % Plotting physical fields
    % plot field, abs
    figure;
    imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi_vs_mode(:,:,ii) ) );
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



























