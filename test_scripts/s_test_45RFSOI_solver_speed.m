% authors: bohan zhang
%
% testing solver speed vs. discretization size

clear; close all;

% dependencies
addpath(genpath('..'));                                                     % all repo codes

% initial settings
disc                = 10;
units               = 'nm';
lambda              = 1300;
index_clad          = 1.0;              % 1.448;
domain              = [2560, 800];      % [y,x]
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

% -------------------------------------------------------------------------
% versus discretization size
% -------------------------------------------------------------------------

% % pick discretizations to try
% discs = [ 80, 40, 20, 10, 5 ];
% 
% % simulation settings
% num_modes   = 10;
% BC          = 0;
% pml_options = [ 1, 200, 20, 2 ];
% guessk      = 0.0098 + 1i * 2.428 * 1e-04;
% 
% % elapsed time
% comp_times_vs_disc = zeros( size(discs) );
% 
% for ii = 1:length(discs)
%     % sweep and time solver
%     
%     fprintf('loop %i of %i\n', ii, length(discs) );
%     
%     % make object
%     Q = c_synthGrating( 'discretization',   discs(ii),       ...
%                         'units',            units,      ...
%                         'lambda',           lambda,     ...
%                         'background_index', index_clad, ...
%                         'domain_size',      domain,     ...
%                         'optimal_angle',    optimal_angle,      ...
%                         'coupling_direction', coupling_direction, ...
%                         'data_directory',   data_dir, ...
%                         'data_filename',    data_filename, ...
%                         'data_notes',       data_notes, ...
%                         'data_mode',        'new', ...
%                         'num_par_workers',  n_workers, ...
%                         'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
%                 );
% 
% 
%     % make grating cell
%     GC = Q.h_makeGratingCell(  Q.convertObjToStruct(), ...
%                                 period, ...
%                                 fill_top, ...
%                                 fill_bot, ...
%                                 offset );
% 
%     % run sim
%     tic;
%     GC = GC.runSimulation( num_modes, BC, pml_options, guessk );        
%     comp_times_vs_disc(ii) = toc
%     
% end
% 
% % plot computational time versus discretization size
% figure;
% plot( discs, comp_times_vs_disc, '-o' );
% xlabel('discretization size'); ylabel('time');
% title('Computational solver time vs. discretization size');
% makeFigureNice();
% 
% % plot slope
% dcomp_time_ddisc = diff( comp_times_vs_disc )./diff( discs );
% figure;
% plot( 1:length(dcomp_time_ddisc), dcomp_time_ddisc, '-o' );
% xlabel('discretization size'); ylabel('slope');
% title('deriv Computational solver time vs. discretization size');
% makeFigureNice();

% -------------------------------------------------------------------------
% versus # of modes
% -------------------------------------------------------------------------

% pick discretizations to try
disc = 10;

% simulation settings
num_modes   = 1:5:51;
BC          = 0;
pml_options = [ 1, 200, 20, 2 ];
guessk      = 0.0098 + 1i * 2.428 * 1e-04;

% elapsed time
comp_times_vs_num_modes = zeros( size(num_modes) );

for ii = 1:length(num_modes)
    % sweep and time solver
    
    fprintf('loop %i of %i\n', ii, length(num_modes) );
    
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

    % run sim
    tic;
    GC = GC.runSimulation( num_modes(ii), BC, pml_options, guessk );        
    comp_times_vs_num_modes(ii) = toc
    
end

% plot computational time versus discretization size
figure;
plot( num_modes, comp_times_vs_num_modes, '-o' );
xlabel('num modes'); ylabel('time');
title('Computational solver time vs. num of modes');
makeFigureNice();

% % plot slope
% dcomp_time_ddisc = diff( comp_times_vs_num_modes )./diff( discs );
% figure;
% plot( 1:length(dcomp_time_ddisc), dcomp_time_ddisc, '-o' );
% xlabel('num modes'); ylabel('slope');
% title('deriv Computational solver time vs. discretization size');
% makeFigureNice();