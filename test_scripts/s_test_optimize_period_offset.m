% authors: bohan zhang
%
% testing 2D optimizer
% for a given set of FF's, find the period and offset that gives desired
% angle and best directivity in two ways:
%   brute force sweep
%   using an optimizer

clear; close all;

% add path to the 45RFSOI functions
addpath(['..' filesep '45RFSOI']);
% add main code
addpath([ '..' filesep 'main' ]);

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
period_vec      = [700, 900];
offset_vec      = linspace(0, 0.3, 2);
fill_top_vec    = 0.5;
fill_bot_vec    = 0.5;

% number of parallel workers, unused
n_workers = 4;

% waveguide index/thickness, unused
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

% chosen fills
fill_top = 0.8;
fill_bot = 0.8;

% period and offset ranges to sweep
periods = 600:10:900;
offsets = 0:0.02:0.98;
% % TEMP
% periods = 00;
% offsets = [0, 0.5];

% simulation settings
num_modes   = 5;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 5, 2 ];

% init saving variables
directivities   = zeros( length(periods), length(offsets) );    % up/down directivity, dimensions period vs. offset
angles_up       = zeros( length(periods), length(offsets) );    % up angle, dimensions period vs. offset
angles_down     = zeros( length(periods), length(offsets) );    % down angle, dimensions period vs. offset

i_loop = 0;
tic;

% run loops
for i_period = 1:length(periods)
    
    for i_offset = 1:length(offsets)
        
        % print loop #
        i_loop = i_loop + 1;
        fprintf('Loop %i of %i\n', i_loop, length(periods)*length(offsets) );
        
        % make grating coupler object
        GC = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), periods(i_period), fill_top, fill_bot, offsets(i_offset) );
        
        % DEBUG plot index
        GC.plotIndex();
        
        % run sim
        GC = GC.runSimulation( num_modes, BC, pml_options );
        
        % DEBUG plot field
        GC.plotEz_w_edges();
       
        % save results
        directivities( i_period, i_offset ) = GC.directivity;
        angles_up( i_period, i_offset )     = GC.max_angle_up;
        angles_down( i_period, i_offset )   = GC.max_angle_down;
        
        toc;
        
    end
    
end
        
% % test the make 45RFSOI function
% period          = 700;
% fill_top        = 0.8;
% fill_bot        = 0.8;
% offset_ratio    = 0.3;
% GC              = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );
% 
% % plot index
% GC.plotIndex();
%         
% % simulate
% % run simulation
% num_modes   = 5;
% BC          = 0;     % 0 for PEC, 1 for PMC
% % PML_options(1): PML in y direction (yes=1 or no=0)
% % PML_options(2): length of PML layer in nm
% % PML_options(3): strength of PML in the complex plane
% % PML_options(4): PML polynomial order (1, 2, 3...)
% pml_options = [ 1, 200, 20, 2 ];
% 
% % run simulation
% tic;
% GC = GC.runSimulation( num_modes, BC, pml_options );
% toc;
% 
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
        
        
        
        
        
        
        
        
        
        
        
        
        