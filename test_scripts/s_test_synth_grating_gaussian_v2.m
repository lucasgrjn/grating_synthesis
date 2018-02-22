% authors: bohan
% 
% script for testing the new, NEW synthesis pipeline object

clear; close all;

% dependencies
addpath(['..' filesep 'main']);                                             % main code
addpath(['..' filesep '45RFSOI']);                                          % 45rf functions

% initial settings
disc                = 10;
units               = 'nm';
lambda              = 1550;
index_clad          = 1.0;
domain              = [2000, 800];      % useful
optimal_angle       = 15;             % still useful
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
                    'num_par_workers',  n_workers ...
            );
%             'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...

        
% run synthesize gaussian grating
% note units are in 'nm'
% 
% angle   = 20; 
MFD     = 10000;
DEBUG   = true;
Q       = Q.synthesizeGaussianGrating(MFD, DEBUG);







% directivity vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, 10*log10(Q.directivities_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Directivity (dB) vs. fill factors');
savefig('dir_v_ff.fig');
saveas(gcf, 'dir_v_ff.png');

% angles vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, Q.angles_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Angles (deg) vs. fill factors');
savefig('angle_v_ff.fig');
saveas(gcf, 'angle_v_ff.png');

% scattering strength alpha vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, real(Q.scatter_str_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Scattering strength (real) vs. fill factors');
savefig('scatter_str_v_ff.fig');
saveas(gcf, 'scatter_str_v_ff.png');

% period vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, Q.periods_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title(['Period (' Q.units.name ') vs. fill factors']);
savefig('period_v_ff.fig');
saveas(gcf, 'period_v_ff.png');

% offset vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, Q.offsets_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Offset vs. fill factors');
savefig('offsets_v_ff.fig');
saveas(gcf, 'offsets_v_ff.png');

% k vs. fill
figure;
imagesc( Q.fill_bots, Q.fill_tops, real(Q.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Real k vs. fill factors');
savefig('k_real_v_ff.fig');
saveas(gcf, 'k_real_v_ff.png');

figure;
imagesc( Q.fill_bots, Q.fill_tops, imag(Q.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Imag k vs. fill factors');
savefig('k_imag_v_ff.fig');
saveas(gcf, 'k_imag_v_ff.png');





        