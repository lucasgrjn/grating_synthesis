% authors: bohan
% 
% script for testing the new synthesis pipeline object
% synth and simulate a uniform grating from previously loaded data

clear; close all;

% % initial settings
% disc        = 10;
% units       = 'nm';
% lambda      = 1550;
% index_clad  = 1.0;
% domain      = [ 1600, 800 ];

% directory to save data to
data_dir        = [ pwd, filesep, 'test_datasave' ];
data_filename   = '2017_11_07 23_19_26 test synth grating.mat';

% make object
Q = c_synthGrating( 'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_mode',        'load' ...
            );


% DEBUG testing synthesis of uniform grating
angle   = 15; 
MFD     = 10000;
Q = Q.synthesizeUniformGrating(angle, MFD);







        
        