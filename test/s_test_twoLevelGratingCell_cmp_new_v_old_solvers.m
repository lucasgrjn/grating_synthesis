% authors: bohan

% Script for testing and debugging the new two level grating cell class
% comparing the new and old solvers

clear; close all;

% dependencies
addpath(['..' filesep 'main']);
addpath(['..' filesep 'auxiliary_functions']);
addpath(['..' filesep 'eme' ]);

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
domain      = [ 3000, 800 ];
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
% wgs_offsets     = [ 0, offset*period ];
wgs_offsets     = [ 0, 200 ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, wg_index, ...
                                        wgs_duty_cycles, wgs_offsets );
                                 
                                 
% DEBUG plot the index
GC.plotIndex();

% -------------------------------------------------------------------------
% Run new solver
% -------------------------------------------------------------------------

% run simulation
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
guessk      = 2*pi*wg_index(1)*fill/lambda;
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 20, 2 ];

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

% DEBUG plot physical fields and all fields
k_all       = GC.debug.k_all;
neff_all    = k_all/(2*pi/lambda);

% % plot all fields
% for ii = 1:length( k_all )
%     % Plotting physical fields
%     % plot field, abs
%     figure;
%     imagesc( GC.x_coords, GC.y_coords, abs( GC.debug.phi_all(:,:,ii) ) );
%     colorbar;
%     set( gca, 'YDir', 'normal' );
%     title( sprintf( 'Field (abs) for mode %i, k = %f + i%f', ii, real( k_all(ii) ), imag( k_all(ii) ) ) );
% end

% plot real and imag k
k_labels = {};
for ii = 1:length( k_all )
    k_labels{end+1} = [ ' ', num2str(ii) ];
end
figure;
plot( real( k_all ), imag( k_all ), 'o' ); 
text( real( k_all ), imag( k_all ), k_labels );
xlabel('real k'); ylabel('imag k');
title('real vs imag k');
makeFigureNice();

% Plot the accepted mode
figure;
imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( ['Field (abs) for accepted mode, k= ', num2str( GC.k  ) ] );

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

% plot all modes
f_plot_all_modes_gui( GC.debug.phi_all, GC.x_coords, GC.y_coords, GC.debug.k_all );



% -------------------------------------------------------------------------
% Run old solver
% -------------------------------------------------------------------------

% run simulation
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
guessk      = 2*pi*wg_index(1)*fill/lambda;
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 500, 2 ];

% run simulation
GC_old = GC.runSimulation_old( num_modes, BC, pml_options, guessk );

% DEBUG plot physical fields and all fields
k_all_old   = GC_old.debug.k_all;
neff_all    = k_all/(2*pi/lambda);

% for ii = 1:length( k_all )
%     % Plotting physical fields
%     % plot field, abs
%     figure;
%     imagesc( GC_old.x_coords, GC_old.y_coords, abs( GC_old.debug.phi_all(:,:,ii) ) );
%     colorbar;
%     set( gca, 'YDir', 'normal' );
%     title( sprintf( 'Field (abs) for mode %i, k = %f + i%f, old solver', ii, real( k_all_old(ii) ), imag( k_all_old(ii) ) ) );
% end

% plot real and imag k
k_labels = {};
for ii = 1:length( k_all_old )
    k_labels{end+1} = [ ' ', num2str(ii) ];
end
figure;
plot( real( k_all_old ), imag( k_all_old ), 'o' ); 
text( real( k_all_old ), imag( k_all_old ), k_labels );
xlabel('real k'); ylabel('imag k');
title('real vs imag k, old solver');
makeFigureNice();

% Plot the accepted mode
figure;
imagesc( GC_old.x_coords, GC_old.y_coords, abs( GC_old.Phi ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( ['Field (abs) for accepted mode, old solver, k = ', num2str( GC_old.k  ) ] );

% display calculated k
fprintf('\nOld Complex k = %f + %fi\n', real(GC_old.k), imag(GC_old.k) );

% display radiated power
fprintf('\nOld Rad power up = %e\n', GC_old.P_rad_up);
fprintf('Old Rad power down = %e\n', GC_old.P_rad_down);
fprintf('Old Up/down power directivity = %f\n', GC_old.directivity);

% display angle of radiation
fprintf('\nOld Angle of maximum radiation = %f deg\n', GC_old.max_angle_up);

% plot full Ez with grating geometry overlaid
GC_old.plotEz_w_edges();
axis equal;

% plot all modes
f_plot_all_modes_gui( GC_old.debug.phi_all, GC_old.x_coords, GC_old.y_coords, GC_old.debug.k_all );


% -------------------------------------------------------------------------
% Run fmm/eme
% -------------------------------------------------------------------------

% Run it in EME
% Set Up Simulation
% note that emeSim uses 'z' as propagation direction and 'x'
% as transverse (synthGrating uses 'x' and 'y' respectively)
% and units are in um
n_cells     = 10;
um          = 1e6;
nm          = 1e9;
dx          = disc*um/nm;                                               % in um
dz          = dx;                                                       % in um
pol         = 0;                                                        % 0 for TE, 1 for TM
xf          = domain(1)*um/nm;                                          % in um
zf          = n_cells * domain(2) * um/nm;                              % in um
lambda_um   = lambda * um/nm;                                           % wl in um
eme_obj     = emeSim(   'discretization', [dx dz], ...
                    'pml', 0.2, ...
                    'domain', [xf zf], ...
                    'backgroundIndex', 1, ...
                    'wavelengthSpectrum', [lambda_um lambda_um 0.1], ...
                    'debug', 'no',...                   
                    'polarization', pol );
diel        = eme_obj.diel;
% grab emeSim coordinates
z_coords_eme    = eme_obj.domain.z;
cur_z           = z_coords_eme(1);          % current z coordinate

% draw grating
eme_obj.diel = repmat( GC.N, 1, n_cells );

% run EME sim
% Converts the dielectric distribution into layers for eigen mode expansion
eme_obj = eme_obj.convertDiel();   
% Runs simulation
eme_obj = eme_obj.runSimulation('plotSource','yes');      
% compute fiber overlap
MFD = 10;
eme_obj = eme_obj.fiberOverlap( 'zOffset', 0:.1:12,...
                                'angleVec', -45:1:45,...
                                'MFD',      MFD,...
                                'overlapDir', 'down' );
% DEBUG show results
gratingUI(eme_obj);














