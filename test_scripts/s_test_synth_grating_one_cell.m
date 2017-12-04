% authors: bohan
% 
% script for testing the new synthesis pipeline object
% tests a single grating cell

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% directory to save data to
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% THIS DOESNT MATTER sweep parameters
period_vec = [700, 900];
offset_vec = linspace(0, 0.3, 2);
ratio_vec  = linspace(0.7, 1.0, 1);
fill_vec   = linspace(0.5, 0.8, 1);

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
                    'ratio_vec',        ratio_vec,  ...
                    'fill_vec',         fill_vec,   ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers ...
            );


% test one grating cell
period  = 780;
fill    = 0.82;
ratio   = 0.7714;
offset  = 0.1163;
[Q, GC] = Q.testMakeGratingCell( period, fill, ratio, offset )

% plot stuff
GC.plotIndex();
GC.plotEz();







































        
        
