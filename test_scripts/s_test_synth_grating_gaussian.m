% authors: bohan
% 
% script for testing the new synthesis pipeline object
% synth and simulate a gaussian grating from prev. simulated data

clear; close all;

% directory to save data to
% data_dir        = [ pwd, filesep, 'test_datasave' ];
% % 11/28/17 data -------------
% data_dir        = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave\2017 11 28 sweep fill\sweep_fill';
% data_filename   = 'sweep fill combined.mat';
% 11/29/17 data -------------
data_dir        = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave\2017 11 29 sweep fill\sweep_fill';
data_filename   = 'sweep fill combined.mat';

% make object
Q = c_synthGrating( 'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_mode',        'load' ...
            );


% run synthesize gaussian grating
% note units are in 'nm'
% 
angle   = 20; 
MFD     = 10000;
Q       = Q.synthesizeGaussianGrating(angle, MFD);






        