% authors: bohan zhang
%
% script for plotting bandstructure vs. period of the 45RFSOI cell

clear; close all;

% add path to the 45RFSOI functions
addpath(['..' filesep '45RFSOI']);
% add main code
addpath([ '..' filesep 'main' ]);

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

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
% unused in this script
period_vec = [700, 900];
offset_vec = linspace(0, 0.3, 2);
% ratio_vec  = linspace(0.7, 1.0, 1);
% fill_vec   = linspace(0.5, 0.8, 1);
fill_top_vec = 0.5;
fill_bot_vec = 0.5;

% number of parallel workers
n_workers = 4;

% waveguide index/thickness
waveguide_index     = [ 3.47, 3.47 ];
waveguide_thicks    = [ 100, 100 ];

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
                    'period_vec',       period_vec, ...
                    'offset_vec',       offset_vec, ...
                    'fill_top_vec',     fill_top_vec, ...
                    'fill_bot_vec',     fill_bot_vec, ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
            );

% init some variables
all_k       = [];
all_periods = [];
chosen_k    = [];
angles      = [];
        
% Loop through different periods and plot the bandstructure vs. period
offset_ratio    = 0.0;
fill_top        = 0.0;
fill_bot        = 0.3;
periods         = 200:20:1500;

tic;
for ii = 1:length(periods)
    
    fprintf('Loop %i of %i\n', ii, length(periods) );
    
    % make grating cell
    GC = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), periods(ii), fill_top, fill_bot, offset_ratio );
    
    % simulate
    % run simulation
    num_modes   = 5;
    BC          = 0;     % 0 for PEC, 1 for PMC
    % PML_options(1): PML in y direction (yes=1 or no=0)
    % PML_options(2): length of PML layer in nm
    % PML_options(3): strength of PML in the complex plane
    % PML_options(4): PML polynomial order (1, 2, 3...)
    pml_options = [ 1, 200, 500, 2 ];

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options );
    k  = GC.k;
    if real(k) <= 0
        % correction for backwards prop mode
        k = -k;
    end
    
    % save all k's
    all_k               = [ all_k, GC.debug.k_all.' ];
    all_periods         = [ all_periods, repmat( periods(ii), 1, length(GC.debug.k_all) ) ];
    chosen_k(end+1)     = k;
    
    % chosne angles
    angles(end+1) = GC.max_angle_down;
    
    toc;
    
end

% plot the bandstructure for ALL modes
figure;
plot( real(all_k), all_periods, 'o' ); hold on;
plot( imag(all_k), all_periods, 'o' ); hold on;
xlabel('k'); ylabel('period (nm)');
legend('real', 'imag');
title('Bandstructure of all solved modes');
makeFigureNice();

% plot the bandstructure
figure;
plot( real(chosen_k), periods, 'o' ); hold on;
plot( imag(chosen_k), periods, 'o' ); hold on;
xlabel('k'); ylabel('period (nm)');
legend('real', 'imag');
title('Bandstructure of chosen modes');
makeFigureNice();

% plot angles
figure;
plot( periods, angles, '-o' );
xlabel('period (nm)'); ylabel('angle (deg)');
title('Angle vs. period');
makeFigureNice();
        
% % test the make 45RFSOI function
% period          = 700;
% fill_top        = 0.8;
% fill_bot        = 0.8;
% offset_ratio    = 0.3;
% GC              = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );
% 
% % plot index
% GC.plotIndex();
%         
% % simulate
% % run simulation
% num_modes   = 5;
% BC          = 0;     % 0 for PEC, 1 for PMC
% % PML_options(1): PML in y direction (yes=1 or no=0)
% % PML_options(2): length of PML layer in nm
% % PML_options(3): strength of PML in the complex plane
% % PML_options(4): PML polynomial order (1, 2, 3...)
% pml_options = [ 1, 200, 500, 2 ];
% 
% % run simulation
% GC = GC.runSimulation( num_modes, BC, pml_options );
% 
% % Plot the accepted mode
% figure;
% imagesc( GC.x_coords, GC.y_coords, abs( GC.Phi ) );
% colorbar;
% set( gca, 'YDir', 'normal' );
% title( sprintf( 'Field (abs) for accepted mode, ka/2pi real = %f', real( GC.k*domain(2)/(2*pi) ) ) );
% 
% % display calculated k
% fprintf('\nComplex k = %f + %fi\n', real(GC.k), imag(GC.k) );
% 
% % display radiated power
% fprintf('\nRad power up = %e\n', GC.P_rad_up);
% fprintf('Rad power down = %e\n', GC.P_rad_down);
% fprintf('Up/down power directivity = %f\n', GC.directivity);
% 
% % display angle of radiation
% fprintf('\nAngle of maximum radiation = %f deg\n', GC.max_angle_up);
% 
% % plot full Ez with grating geometry overlaid
% GC.plotEz_w_edges();
% axis equal;
        
        
        
        
        
        
        
        
        
        
        
        
        