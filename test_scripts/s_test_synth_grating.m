% authors: bohan
% 
% script for testing the new synthesis pipeline object

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 1600, 800 ];

% directory to save data to
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
% % case 1
% period_vec = [ 800, 850 ];
% offset_vec = [ 0, 0.2 ];
% ratio_vec  = [ 1, 0.8 ];
% fill_vec   = [ 0.75, 0.75 ];
% case 2
% period_vec = [ 1, 2 ];
% offset_vec = [ 3, 4 ];
% ratio_vec  = [ 5, 6];
% fill_vec   = [ 7, 8 ];
% % case 3 only change period
% period_vec = 600:100:1000;
% offset_vec = 0;
% ratio_vec  = 1;
% fill_vec   = 0.5;
% % case 4 only change offset
% period_vec = 600;
% offset_vec = 0:100:400;
% ratio_vec  = 1;
% fill_vec   = 0.5;
% % case 5 only change ratio
% period_vec = 1000;
% offset_vec = 0;
% ratio_vec  = linspace(0.5, 1, 24*10);
% fill_vec   = 0.5;
% % case 6 only change fill
% period_vec = 600;
% offset_vec = 0;
% ratio_vec  = 1;
% fill_vec   = linspace(0.5, 1, 2);
% % case 7 one case
% period_vec = 600;
% offset_vec = 0;
% ratio_vec  = 1;
% fill_vec   = 0.7;
% case 8 a mini sweep
period_vec = 1000:200:1200;
offset_vec = linspace(0, 0.3, 2);
ratio_vec  = linspace(0.7, 1.0, 1);
fill_vec   = linspace(0.5, 0.8, 1);
% % case 9 a larger sweep
% period_vec = linspace(500, 800, 13);
% offset_vec = linspace(0, 0.9, 10);
% ratio_vec  = linspace(0.1, 1.1, 11);
% fill_vec   = linspace(0.1, 0.9, 9);
% % case 10 an even larger sweep
% period_vec = linspace(500, 800, 15);
% offset_vec = linspace(0, 0.9, 15);
% ratio_vec  = linspace(0.1, 1.3, 15);
% fill_vec   = linspace(0.1, 1.0, 15);
% % case 11 a larger sweep with evenly spaced periods
% period_vec = linspace(500, 800, 16);
% offset_vec = linspace(0, 0.9, 10);
% ratio_vec  = linspace(0.1, 1.1, 11);
% fill_vec   = linspace(0.1, 0.9, 9);
% % case 12 even larger sweep
% period_vec = 500:10:800;
% offset_vec = linspace(0, 0.95, 20);
% ratio_vec  = linspace(0.1, 1.2, 20);
% fill_vec   = linspace(0.1, 0.95, 20);

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
        
% run parameter sweep
% Q = Q.runParameterSweep();

% DEBUG make a single grating cell
% Q = Q.testMakeGratingCell( period_vec(2), fill_vec, ratio_vec, offset_vec(1) );
% Q = Q.testMakeGratingCell( period_vec(2), fill_vec, ratio_vec, 0.4 );

% DEBUG running the fill/period angle sweep
Q = Q.sweepPeriodFill();









































        
        
