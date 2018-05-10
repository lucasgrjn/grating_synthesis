% authors: bohan zhang
%
% testing the two solver versions, one non symmetric the other is hermitian
% (without pmls)

clear; close all;

% dependencies
addpath(genpath('..'));                                                     % all repo codes

% initial settings
disc                = 10;
units               = 'nm';
lambda              = 1300;
index_clad          = 1.0;              % 1.448;
domain              = [2500, 800];      % [y,x]
optimal_angle       = 20;               
coupling_direction  = 'down';
% settings below are not currently used
data_dir            = 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\data';
data_filename       = 'lol.mat';
data_notes          = 'meh';
data_mode           = 'new';
n_workers           = 1;

% grating geometry
period              = 560;
fill_top_bot_ratio  = 0.325;
fill_bot            = 0.6;
fill_top            = fill_bot * fill_top_bot_ratio;
offset              = 0.9967;

% simulation settings
num_modes   = 10;
BC          = 0;
pml_options = [ 1, 200, 20, 2 ];
guessk      = 0.0098 + 1i * 2.428 * 1e-04;


% make object
Q = c_synthGrating( 'discretization',   disc,       ... 
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'optimal_angle',    optimal_angle,      ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
            );


% make grating cell
GC = Q.h_makeGratingCell(  Q.convertObjToStruct(), ...
                            period, ...
                            fill_top, ...
                            fill_bot, ...
                            offset );

% run sim with new hermitian modesolver
tic;
GC2 = GC.runSimulation_v2_symm( num_modes, BC, pml_options, guessk );        
toc;
                        
% run sim with original non-symmetric modesolver
tic;
GC1 = GC.runSimulation( num_modes, BC, pml_options, guessk );        
toc;

GC1.plot_E_field_gui();
GC2.plot_E_field_gui();



































