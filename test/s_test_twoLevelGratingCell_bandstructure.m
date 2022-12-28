% authors: bohan
% Calculates bandstructure for two level grating

clear; close all;

% dependencies
addpath(['..' filesep 'main']);
addpath(['..' filesep '45RFSOI']);

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2500, 800 ];

% directory to save data to
% unused for this script
data_dir        = '';
data_filename   = '';
data_notes      = '';

% number of parallel workers, unused
n_workers = 0;

% desired angle
optimal_angle = 15;

% coupling up/down
coupling_direction = 'down';

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

      


% -------------------------------------------------------------------------
% Calc bandstructure
% -------------------------------------------------------------------------


% grating coupler object dimensions
period          = 690;
offset_ratio    = 0.2;
fill_bot        = 0.9;
fill_top        = fill_bot * 0.9;

% simulation settings
num_modes   = 1;
BC          = 0;                                                                % 0 = PEC
pml_options = [1, 200, 20, 2]; 
guessk      = 1i * 5.404e-5 + 0.009631;
OPTS        = struct( 'fix_neg_k', true );

% % make grating coupler object
% GC = Q.h_makeGratingCell( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );
% 
% % run sim
% GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
% % display stuff
% GC.plotIndex();
% GC.plotEz_w_edges();


% Pick frequencies/wavelengths to sweep
lambda_max          = lambda * 1.2;
lambda_min          = lambda / 1.2;
n_lambda            = 30;
lambda_all_smaller  = linspace(lambda_min, lambda, n_lambda/2);
lambda_all_larger   = linspace(lambda, lambda_max, n_lambda/2);
lambda_all          = [ lambda_all_smaller(1:end-1), lambda_all_larger ];
% k0_max      = 10/period;
% k0_min      = 0.5/period;
% k0_all      = linspace( k0_min, k0_max, 200 );
% lambda      = 2*pi./k0_all;
% lambda_min  = 2*pi/k0_max;
% lambda_max  = 2*pi/k0_min;
% lambda      = linspace( lambda_min, lambda_max, 100 );

% init saving vars
% k_all               = zeros( size(lambda_all) );
% angle_down          = zeros( size(lambda_all) );
k_all_larger        = zeros( size(lambda_all_larger) );
angle_down_larger   = zeros( size(lambda_all_larger) );
k_all_smaller       = zeros( size(lambda_all_smaller) );
angle_down_smaller  = zeros( size(lambda_all_smaller) );

% Calculate bandstructure in two steps

% first sweep larger lambda
tic;
for ii = 1:length(lambda_all_larger)
    
    fprintf('loop %i of %i\n', ii, length(lambda_all_larger) );
    
    % set wavelength
    Q.lambda = lambda_all_larger(ii);
    
    % make grating coupler object
    GC = Q.h_makeGratingCell( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );

    % run sim
    GC = GC.runSimulation( num_modes, BC, pml_options, guessk, OPTS );
    
    % save stuff
    k_all_larger(ii)       = GC.k;
    angle_down_larger(ii)  = GC.max_angle_down;
    
    % update guessk
    guessk = GC.k;
    
    toc;
    
end

% then sweep smaller lambda in backwards order
guessk = k_all_larger(1);
for ii = length(lambda_all_smaller):-1:1
    
    fprintf('loop %i of %i\n', ii, length(lambda_all_smaller) );
    
    % set wavelength
    Q.lambda = lambda_all_smaller(ii);
    
    % make grating coupler object
    GC = Q.h_makeGratingCell( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );

    % run sim
    GC = GC.runSimulation( num_modes, BC, pml_options, guessk, OPTS );
    
    % save stuff
    k_all_smaller(ii)       = GC.k;
    angle_down_smaller(ii)  = GC.max_angle_down;
    
    % update guessk
    guessk = GC.k;
    
    toc;
    
end


% finally stitch both sweeps together
k_all       = [ k_all_smaller(1:end-1), k_all_larger ];
angle_down  = [ angle_down_smaller(1:end-1), angle_down_larger ];


% plot bandstructure, real k vs. lambda
figure;
plot( lambda_all, real(k_all), '-o' ); 
xlabel('lambda (nm)'); ylabel('k');
title('real k vs. lambda');
makeFigureNice();

% plot bandstructure, imag k vs. lambda
figure;
plot( lambda_all, imag(k_all), '-o' ); 
xlabel('lambda (nm)'); ylabel('k');
title('imag k vs. lambda');
makeFigureNice();

% plot angle vs. lambda
figure;
plot( lambda_all, angle_down, '-o' );
xlabel('lambda (nm)'); ylabel('degrees');
title('angle downwards vs. lambda');
makeFigureNice();



% figure;
% plot( real(k_all)*period, (2*pi./lambda)*period, '-o' ); hold on;
% plot( imag(k_all)*period, (2*pi./lambda)*period, '-o' );
% xlabel('ka'); ylabel('k_0a');
% title('band structure');
% makeFigureNice();

% % plot all modes
% figure;
% plot( real(k_all_all)*period, k0_a_allall, 'o' );
% xlabel('ka'); ylabel('k_0a');
% title('band structure of ALL solved modes');
% makeFigureNice();
% 
% % plot bandstructure vs. period @ 1550nm
% % LOL can i do this
% a_all = (2*pi./lambda)*period/( 2*pi/1550 );
% 
% % plot bandstructure
% figure;
% plot( real(k_all)*period, a_all, '-o' ); hold on;
% plot( imag(k_all)*period, a_all, '-o' );
% xlabel('ka'); ylabel('a (period, nm)');
% title('band structure vs. period @ 1550nm');
% makeFigureNice();























