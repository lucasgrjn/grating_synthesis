% authors: bohan zhang
%
% script for testing the new f_makeGratingCell_45RFSOI function

clear; close all;

% dependencies
addpath(['..' filesep '45RFSOI']);                                          % 45rf
addpath([ '..' filesep 'main' ]);                                           % main
addpath([ '..' filesep 'auxiliary_functions' ]);                            % gui

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% directory to save data to
% unused for this script
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
% unused in this script
period_vec = [700, 900];
offset_vec = linspace(0, 0.3, 2);
% ratio_vec  = linspace(0.7, 1.0, 1);
% fill_vec   = linspace(0.5, 0.8, 1);
fill_top_vec = 0.5;
fill_bot_vec = 0.5;

% number of parallel workers
n_workers = 4;

% waveguide index/thickness
waveguide_index     = [ 3.47, 3.47 ];
waveguide_thicks    = [ 100, 100 ];

% desired angle
optimal_angle = 15;

% coupling up/down
coupling_direction = 'down';

% make object
Q = c_synthGrating( 'discretization',   disc,       ...
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'period_vec',       period_vec, ...
                    'offset_vec',       offset_vec, ...
                    'fill_top_vec',     fill_top_vec, ...
                    'fill_bot_vec',     fill_bot_vec, ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
            );
        
% test the make 45RFSOI function
period          = 700;
fill_top        = 0.3;
fill_bot        = 0.3;
offset_ratio    = 0.3;
GC              = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );

% plot index
GC.plotIndex();
        
% simulate
% run simulation
num_modes   = 5;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 20, 2 ];

% guessk
guessk = 2*2*pi/lambda;
% guessk  = 0.0094 + 0.0001i;
% guessk = (0.002634 - 0.000038i);
% guessk = -0.0019 + 0.0000i;

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
fprintf('\nAngle of maximum radiation up = %f deg\n', GC.max_angle_up);
fprintf('\nAngle of maximum radiation down = %f deg\n', GC.max_angle_down);

% plot full Ez with grating geometry overlaid
GC.plotEz_w_edges();
axis equal;
        
% plot all modes
f_plot_all_modes_gui( GC.debug.phi_all, GC.x_coords, GC.y_coords, GC.debug.k_all )
        
        
        
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
zf          = n_cells * period * um/nm;                                 % in um
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
        
        
        
        
        
        
        