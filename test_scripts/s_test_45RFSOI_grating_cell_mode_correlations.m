% authors: bohan zhang
%
% script for testing the new f_makeGratingCell_45RFSOI function
%   testing correlations between modes

clear; close all;

% dependencies
addpath(genpath('..'));                                                     % all repo codes


% initial settings
disc                = 10;
units               = 'nm';
lambda              = 1200;
index_clad          = 1.0; % 1.448;
domain              = [2500, 800];      % useful
optimal_angle       = 20;             % still useful
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
                    'optimal_angle',    optimal_angle,      ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
            );

% grating geometry
period              = 690;
fill_top_bot_ratio  = 0.0;
fill_bot            = 0.7;
fill_top            = fill_bot * fill_top_bot_ratio;
offset              = 0.64;
        
        
% make grating cell
GC = Q.h_makeGratingCell(  Q.convertObjToStruct(), ...
                            period, ...
                            fill_top, ...
                            fill_bot, ...
                            offset );
                        
% simulation settings
num_modes   = 10;
BC          = 0;
pml_options = [ 1, 200, 20, 2 ];
guessk      = 0.01084 + 1i * 0.0001073;

% run sim
tic;
GC = GC.runSimulation( num_modes, BC, pml_options, guessk );        
toc;

% plot index
GC.plotIndex();

       



% correlate first mode with the other modes
% first_mode              = GC.Phi_vs_mode(:,:,1) .* repmat( exp( 1i * GC.x_coords * real(GC.k_vs_mode( 1 )) ), size(GC.Phi_vs_mode(:,:,1), 1), 1 );
second_mode             = GC.Phi_vs_mode(:,:,2) .* repmat( exp( 1i * GC.x_coords * real(GC.k_vs_mode( 2 )) ), size(GC.Phi_vs_mode(:,:,2), 1), 1 );
[GC, max_overlaps]      = GC.calc_mode_overlaps( second_mode );
        
% plot all modes
GC = GC.plot_E_field_gui();

        

        
        
        
        
        
        
        