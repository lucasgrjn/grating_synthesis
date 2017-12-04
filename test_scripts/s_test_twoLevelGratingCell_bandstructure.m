% authors: bohan

% Script for testing and debugging the new two level grating cell class
% Calculates bandstructure for two level grating to extrapolate ideal range
% of periods to sweep

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% directory to save data to
% unused for this script
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% sweep parameters
% unused in this script
period_vec = [700, 900];
offset_vec = linspace(0, 0.3, 2);
ratio_vec  = linspace(0.7, 1.0, 1);
fill_vec   = linspace(0.5, 0.8, 1);

% number of parallel workers
n_workers = 4;

% waveguide index/thickness
waveguide_index     = [ 3.47, 3.47 ];
waveguide_thicks    = [ 100, 100 ];

% desired angle
optimal_angle = 15;

% coupling up/down
coupling_direction = 'down';

% make synth grating object
Q = c_synthGrating( 'discretization',   disc,       ...
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'period_vec',       period_vec, ...
                    'offset_vec',       offset_vec, ...
                    'ratio_vec',        ratio_vec,  ...
                    'fill_vec',         fill_vec,   ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers ...
            );


        


% -------------------------------------------------------------------------
% Calc bandstructure
% -------------------------------------------------------------------------

% grating cell params
period  = 800;
ratio   = 1;
offset  = 0.0;
fill    = 0.40;

% Pick wavelengths to sweep
k0_max      = 10/period;
k0_min      = 0.5/period;
k0_all      = linspace( k0_min, k0_max, 200 );
lambda      = 2*pi./k0_all;
% lambda_min  = 2*pi/k0_max;
% lambda_max  = 2*pi/k0_min;
% lambda      = linspace( lambda_min, lambda_max, 100 );

% k
k_all       = zeros( size(lambda) );
k_all_all   = [];
k0_a        = (2*pi./lambda)*period;
k0_a_allall = [];

% Calculate bandstructure

% start parpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if ~isempty(poolobj)
    % shut down previously made parallel pool
    delete(gcp('nocreate'));
end
parpool('local', 4);

tic;
parfor ii = 1:length(lambda)
    
    fprintf('loop %i of %i\n', ii, length(lambda) );
    
    % make synth grating object
    Q = c_synthGrating( 'discretization',   disc,       ...
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'period_vec',       period_vec, ...
                    'offset_vec',       offset_vec, ...
                    'ratio_vec',        ratio_vec,  ...
                    'fill_vec',         fill_vec,   ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers ...
            );
    
    % set lambda
    Q.lambda = lambda(ii);
    
    % make/run grating cell
    [Q, GC] = Q.testMakeGratingCell( period, fill, ratio, offset );
    
    % flip sign of k if needed
    if real(GC.k) < 0
        GC.k = -GC.k;
    end
    
    % save k
    k_all(ii)   = GC.k;
    k_all_all   = [ k_all_all, GC.debug.k_all.' ];
    k0_a_allall = [ k0_a_allall, (2*pi/lambda(ii))*period*ones(1, length(GC.debug.k_all) ) ];
    
end
toc;

% plot bandstructure
figure;
plot( real(k_all)*period, (2*pi./lambda)*period, '-o' ); hold on;
plot( imag(k_all)*period, (2*pi./lambda)*period, '-o' );
xlabel('ka'); ylabel('k_0a');
title('band structure');
makeFigureNice();

% plot all modes
figure;
plot( real(k_all_all)*period, k0_a_allall, 'o' );
xlabel('ka'); ylabel('k_0a');
title('band structure of ALL solved modes');
makeFigureNice();

% plot bandstructure vs. period @ 1550nm
% LOL can i do this
a_all = (2*pi./lambda)*period/( 2*pi/1550 );

% plot bandstructure
figure;
plot( real(k_all)*period, a_all, '-o' ); hold on;
plot( imag(k_all)*period, a_all, '-o' );
xlabel('ka'); ylabel('a (period, nm)');
title('band structure vs. period @ 1550nm');
makeFigureNice();























