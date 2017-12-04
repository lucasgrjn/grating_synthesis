% authors: bohan
% 
% script for testing my gaussian generation function

clear; close all;

% % initial settings
% disc        = 10;
% units       = 'nm';
% lambda      = 1550;
% index_clad  = 1.0;
% domain      = [ 1600, 800 ];

% directory to load data from
data_dir        = [ pwd, filesep, 'test_datasave' ];
data_filename   = '2017_11_07 23_19_26 test synth grating.mat';

% make object
Q = c_synthGrating( 'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_mode',        'load' ...
            );
        
% testing the gaussian mode generating function
w0      = 5000;
zvec    = 0;
% yvec    = linspace( 0, 1000, 100 );
xvec    = linspace( -w0*4, w0*4, 100 );
theta   = 45;
d0      = 0;
nclad   = 1;
[Q, F]  = Q.fiberModeGaussian(w0, zvec, xvec, theta, d0, nclad);
% 
% % plot Ez
% figure;
% imagesc( xvec, yvec, abs(F.Ez) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ez');
% 
% % plot Ex
% figure;
% imagesc( xvec, yvec, abs(F.Ex) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ex');
% 
% % plot Ey
% figure;
% imagesc( xvec, yvec, abs(F.Ey) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ey');
%         
% 
% 
% % Compare with Cale's code
% lambda_um   = lambda*1e-3;
% k0          = 2*pi/lambda_um;
% w0_um       = w0*1e-3;
% pol         = 0;                % 0 for TE
% yvec_um     = yvec*1e-3;
% xvec_um     = xvec*1e-3;
% d0_um       = d0*1e-3;
% F_c         = fiberModeGaussian(w0_um, k0, pol, yvec_um, xvec_um, theta, d0_um);
% 
% % plot Ez
% figure;
% imagesc( xvec, yvec, abs(F_c.Ez) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ez, Cale''s');
% 
% % plot Ex
% figure;
% imagesc( xvec, yvec, abs(F_c.Ex) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ex, Cale''s');
% 
% % plot Ey
% figure;
% imagesc( xvec, yvec, abs(F_c.Ey) );
% set( gca, 'ydir', 'normal' );
% colorbar;
% title('DEBUG plot of Gaussian Ey, Cale''s');






        
        