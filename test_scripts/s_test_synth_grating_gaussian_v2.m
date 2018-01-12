% authors: bohan
% 
% script for testing the new, NEW synthesis pipeline object

clear; close all;

% dependencies
addpath(['..' filesep 'main']);                                             % main code
addpath(['..' filesep '45RFSOI']);                                          % 45rf functions

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [2000, 800];      % deprecated?
period_vec  = 0;                % deprecated.
offset_vec  = 0;                % deprecated.
fill_top_vec = 0;               % deprecated.
fill_bot_vec = 0;               % deprecated.
optimal_angle = 15;             % still useful
waveguide_index     = 1;           % should be deprecated
waveguide_thicks    = 1;           % should be deprecated
coupling_direction  = 'down';
data_dir            = 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\data';
data_filename       = 'lol.mat';
data_notes          = 'meh';
data_mode           = 'new';
n_workers           = 1;

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

        
% run synthesize gaussian grating
% note units are in 'nm'
% 
angle   = 15; 
MFD     = 10000;
Q       = Q.synthesizeGaussianGrating(angle, MFD);






        