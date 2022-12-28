% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% debugging the PML
% comparing two types of PMLs, polynomial and polynomial strength + sin^2

clear; close all;

% -------------------------------------------------------------------------
% Debug using 2 level grating cell
% -------------------------------------------------------------------------

% Script for testing and debugging the new two level grating cell class

clear; close all;

% import code
addpath(['..' filesep 'main']);         % main
addpath(['..' filesep '45RFSOI']);      % 45rfsoi

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
% domain      = [ 60, 50 ];
domain      = [ 2000, 700 ];
numcells    = 10;

% unused for this script
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% make the directory to save data to, if not already in existence
% unused for this script
mkdir( data_dir );

% sweep parameters
% unused in this script
period_vec      = [700, 900];
offset_vec      = linspace(0, 0.3, 2);
fill_top_vec    = 0.5;
fill_bot_vec    = 0.5;

% number of parallel workers
% unused for this script
n_workers = 4;

% waveguide index/thickness
waveguide_index     = [ 3.47, 3.47 ];
waveguide_thicks    = [ 100, 100 ];

% desired angle
% unused for this script
optimal_angle = 15;

% coupling up/down
% unused for this script
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
        
% test the make 45RFSOI function
period          = domain(2);
fill_top        = 0.75;
fill_bot        = 1.0;
offset_ratio    = 0.0;
GC              = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );
                                 
                                 
% DEBUG plot the index
GC.plotIndex();

% compare bandstructures
% Pick wavelengths to sweep
k0_max      = 10/period;
k0_min      = 0.5/period;
k0_all      = linspace( k0_min, k0_max, 200 );
lambda      = 2*pi./k0_all;

% init saving variables
all_k1      = [];   % returned k from pml type 1
all_k2      = [];   % returned k from pml type 2
all_k0      = [];
all_k_old   = [];
all_k0_old  = [];

% set simulation options
num_modes   = 1;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 5, 2 ];
DEBUG       = false;

% calc bandstructure
guessk1 = 1e-5;
guessk2 = 1e-5;
guessk_old = 1e-5;
for ii = 1:length(k0_all)
    
    fprintf('\nloop %i of %i\n\n', ii, length(k0_all));
    
    % run simulation
    % run new
    % Phi_all has dimensions ny vs. nx vs. mode #
    fprintf('running new solver, pml type 1\n');
    tic;
    [Phi_all1, k_all1, A, B] = complexk_mode_solver_2D_PML( GC.N, ...
                                                               disc, ...
                                                               k0_all(ii), ...
                                                               num_modes, ...
                                                               guessk1, ...
                                                               BC, ...
                                                               [pml_options, 1], ...
                                                               DEBUG );
    toc;
    
    % run new
    % Phi_all has dimensions ny vs. nx vs. mode #
    fprintf('running new solver, pml type 2\n');
    tic;
    [Phi_all2, k_all2, A, B] = complexk_mode_solver_2D_PML( GC.N, ...
                                                               disc, ...
                                                               k0_all(ii), ...
                                                               num_modes, ...
                                                               guessk2, ...
                                                               BC, ...
                                                               [pml_options, 2], ...
                                                               DEBUG );
    toc;
    
    % run old
    fprintf('running old solver\n');
    tic;
    [Phi_1D_old, k_old, Phi_all_old, k_all_old, A_old, B_old] = complexk_mode_solver_2D_PML_old( GC.N, ...
                                                                                   disc, ...
                                                                                   k0_all(ii), ...
                                                                                   num_modes, ...
                                                                                   guessk_old, ...
                                                                                   BC, ...
                                                                                   pml_options );
    toc;
                 
    all_k1(end+1)           = k_all1;       % returned k from pml type 1
    all_k2(end+1)           = k_all2;       % returned k from pml type 2
    all_k_old(end+1)        = k_all_old;    % old k
    
    % set new guessk's
    guessk1         = all_k1(end);
    guessk2         = all_k2(end);
    guessk_old      = all_k_old(end);
    
end


% plot the bandstructure for ALL modes, new solver pml type 1 and 2
figure;
plot( real(all_k1)*period/pi, k0_all*period/pi, 'o' ); hold on;
plot( imag(all_k1)*period/pi, k0_all*period/pi, 'o' );
plot( real(all_k2)*period/pi, k0_all*period/pi, 'o' );
plot( imag(all_k2)*period/pi, k0_all*period/pi, 'o' );
plot( real(all_k_old)*period/pi, k0_all*period/pi, 'o' );
plot( imag(all_k_old)*period/pi, k0_all*period/pi, 'o' );
xlabel('ka/pi'); ylabel('k0*a/pi');
% legend('real', 'imag', 'center of bandgap');
legend('real, pm1', 'imag, pml1', 'real, pml2', 'imag, pml2', 'real, old', 'imag, old');
title('Bandstructure of all solved modes, new ver. pml 1 vs pml 2 vs. old version');
makeFigureNice();





























