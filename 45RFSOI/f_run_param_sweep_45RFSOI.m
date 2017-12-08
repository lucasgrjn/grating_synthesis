function [] = f_run_param_sweep_45RFSOI( n_workers, fill_top )
% authors: bohan
% 
% function for running the parameter sweep
% FOR THE 45RFSOI dimensions
%
% inputs:
% 	TEMPORARY: pick ONE fill to try

% clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% convert fills to strings for naming
fill_str = strrep( num2str(fill_top), '.', 'd' );
% fill2_str = strrep( num2str(fill2), '.', 'd' );

% directory to save data to
data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = [ 'sweep fill top ' fill_str ];
% data_notes      = [ 'sweep fill 1 ' num2str(fill1) ' fill 2 ' num2str(fill2) ];
data_notes      = [ 'sweep fill top ' num2str(fill_top) ];

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
period_vec      = 600:20:1400;
offset_vec      = linspace(0, 1.0, 40);
% fill_bot_vec    = linspace(0.3, 1.0, 40);
fill_bot_vec    = 0.3:0.02:1.0;
% ratio_vec   = linspace(0.3, 1.2, 40);
fill_top_vec    = fill_top;                         % [ fill1, fill2 ];

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
                    'fill_bot_vec',     fill_bot_vec, ...
                    'fill_top_vec',     fill_top_vec, ...
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
%                     'ratio_vec',        ratio_vec,  ...
%                     'fill_vec',         fill_vec,   ...
        
% run parameter sweep
Q = Q.runParameterSweep();


end % end function






        
        
