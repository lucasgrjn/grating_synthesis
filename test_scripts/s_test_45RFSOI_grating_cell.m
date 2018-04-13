% authors: bohan zhang
%
% script for testing the new f_makeGratingCell_45RFSOI function

clear; close all;

% dependencies
addpath(genpath('..'));                                                     % all repo codes
% addpath(['..' filesep '45RFSOI']);                                          % 45rf
% addpath([ '..' filesep 'main' ]);                                           % main
% addpath([ '..' filesep 'auxiliary_functions' ]);                            % gui

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
period              = 630;
fill_top_bot_ratio  = 0.35;
fill_bot            = 0.8;
fill_top            = fill_bot * fill_top_bot_ratio;
offset              = 0.92;
        
        
% make grating cell
GC = Q.h_makeGratingCell(  Q.convertObjToStruct(), ...
                            period, ...
                            fill_top, ...
                            fill_bot, ...
                            offset );
                        
% simualtion settings
num_modes   = 5;
BC          = 0;
pml_options = [ 1, 200, 20, 2 ];
guessk      = 0.01183 + 1i * 1.58e-4;

% run sim
tic;
GC = GC.runSimulation( num_modes, BC, pml_options, guessk );        
toc;

% plot index
GC.plotIndex();


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
%         
% % plot all modes
% f_plot_all_modes_gui( GC.debug.phi_all, GC.x_coords, GC.y_coords, GC.debug.k_all )
        
        
        
% % -------------------------------------------------------------------------
% % Run fmm/eme
% % -------------------------------------------------------------------------
% 
% % Run it in EME
% % Set Up Simulation
% % note that emeSim uses 'z' as propagation direction and 'x'
% % as transverse (synthGrating uses 'x' and 'y' respectively)
% % and units are in um
% n_cells     = 10;
% um          = 1e6;
% nm          = 1e9;
% dx          = disc*um/nm;                                               % in um
% dz          = dx;                                                       % in um
% pol         = 0;                                                        % 0 for TE, 1 for TM
% xf          = domain(1)*um/nm;                                          % in um
% zf          = n_cells * period * um/nm;                                 % in um
% lambda_um   = lambda * um/nm;                                           % wl in um
% eme_obj     = emeSim(   'discretization', [dx dz], ...
%                     'pml', 0.2, ...
%                     'domain', [xf zf], ...
%                     'backgroundIndex', 1, ...
%                     'wavelengthSpectrum', [lambda_um lambda_um 0.1], ...
%                     'debug', 'no',...                   
%                     'polarization', pol );
% diel        = eme_obj.diel;
% % grab emeSim coordinates
% z_coords_eme    = eme_obj.domain.z;
% cur_z           = z_coords_eme(1);          % current z coordinate
% 
% % draw grating
% eme_obj.diel = repmat( GC.N, 1, n_cells );
% 
% % run EME sim
% % Converts the dielectric distribution into layers for eigen mode expansion
% eme_obj = eme_obj.convertDiel();   
% % Runs simulation
% eme_obj = eme_obj.runSimulation('plotSource','yes');      
% % compute fiber overlap
% MFD = 10;
% eme_obj = eme_obj.fiberOverlap( 'zOffset', 0:.1:12,...
%                                 'angleVec', -45:1:45,...
%                                 'MFD',      MFD,...
%                                 'overlapDir', 'down' );
% % DEBUG show results
% gratingUI(eme_obj);
        
        
        
        
        
        
        