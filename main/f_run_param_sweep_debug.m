function [] = f_run_param_sweep_debug( n_workers, fill )
% authors: bohan
% 
% function for running the parameter sweep
%
% inputs:
% 	TEMPORARY: pick ONE fill to try
%  this one is a small sweep for debugging

% clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% convert fills to strings for naming
fill_str = strrep( num2str(fill), '.', 'd' );

% directory to save data to
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_dir        = 'C:\Users\beezy\Google Drive\research\popovic group\projects\grating synthesis\data\test_datasave';
data_filename   = [ 'sweep fill ' fill_str ];
data_notes      = [ 'sweep fill ' num2str(fill) ];

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
period_vec  = 800:10:840;
offset_vec  = 0;
% ratio_vec   = 1.0;
% fill_vec    = fill;                         % [ fill1, fill2 ];
fill_top_vec = 0.6;
fill_bot_vec = 0.8;

% for debugging
% period_vec = 500:10:800;
% offset_vec = 0.2;
% ratio_vec  = 1.0;
% fill_vec   = 0.6;
% period_vec  = 500;
% offset_vec  = 0.2;
% ratio_vec   = 1.0;


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
                    'num_par_workers',  n_workers ...
            );
%                     'ratio_vec',        ratio_vec,  ...
%                     'fill_vec',         fill_vec,   ...
        
% run parameter sweep
Q = Q.runParameterSweep();


end % end function






        
        
