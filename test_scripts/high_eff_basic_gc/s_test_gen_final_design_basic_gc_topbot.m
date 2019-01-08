% authors: bohan
% 
% script for testing the latest synthesis pipeline
% generate final design using the top/bot sweep

clear; close all;

% dependencies
% % on lab desktop
addpath(genpath('C:\Users\bz\git\grating_synthesis'));                      % grating synthesis codes
% % on laptop
addpath(genpath('C:\Users\beezy\git\grating_synthesis'));                   % grating synthesis codes

% synthesis object to load
% 100nm thick layers, -20 deg, initial test run
filename = 'synth_obj_2019_01_06_12_11_22_lambda1310_optangle-20_NO_GC.mat';
filepath = 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\data\2019 01 19 basic grating topbot\2019_01_06_12_11_22_lambda1310_optangle-20_thick100';

% load synth_obj
load( [ filepath filesep filename ] );

% synthesize final design
MFD             = 10.4 * 1e3;                                                         % in nm
input_wg_type   = 'bottom';
synth_obj       = synth_obj.generate_final_design_gaussian_topbot( MFD, input_wg_type );

% plot final design

% directivity vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.dir), log10(synth_obj.synthesized_design.dir), '-o' );
xlabel('unit cell #'); ylabel('final directivity (log10)');
title('Directivity vs unit cell');
makeFigureNice();

% bottom fill vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.bot_fill), synth_obj.synthesized_design.bot_fill, '-o' );
xlabel('unit cell #'); ylabel('bottom fill');
title('Bottom fill vs unit cell');
makeFigureNice();

% top fill vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.top_fill), synth_obj.synthesized_design.top_fill, '-o' );
xlabel('unit cell #'); ylabel('top fill');
title('Top fill vs unit cell');
makeFigureNice();

% top/bot fill ratio vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.top_bot_fill_ratio), synth_obj.synthesized_design.top_bot_fill_ratio, '-o' );
xlabel('unit cell #'); ylabel('top/bot fill ratio');
title('Top/bot fill ratio vs unit cell');
makeFigureNice();

% period vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '-o' );
xlabel('unit cell #'); ylabel('period');
title('Period vs unit cell');
makeFigureNice();

% offset vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.offset), synth_obj.synthesized_design.offset, '-o' );
xlabel('unit cell #'); ylabel('offset');
title('Offset vs unit cell');
makeFigureNice();

% angle vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '-o' );
xlabel('unit cell #'); ylabel('angle (deg)');
title('Angle vs unit cell');
makeFigureNice();

% scattering strength vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.des_scatter), synth_obj.synthesized_design.des_scatter, '-o' ); hold on;
plot( 1:length(synth_obj.synthesized_design.scatter_str), synth_obj.synthesized_design.scatter_str, '-o' );
xlabel('unit cell #'); ylabel('scattering strength');
legend('Desired','Synthesized');
title('Scattering strength vs unit cell');
makeFigureNice();

% final index distribution
figure;
imagesc( synth_obj.synthesized_design.N );
set( gca, 'ydir', 'normal' );
colorbar;
title('Final index distribution');

% save final object
save( [ filepath filesep filename(1:end-4) '_final.mat' ], 'synth_obj' );















