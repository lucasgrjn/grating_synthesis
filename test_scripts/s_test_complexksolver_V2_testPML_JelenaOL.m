% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% debugging the PML
% comparing two types of PMLs, polynomial and polynomial strength + sin^2
% recreating the results from jelena's OL paper (fig 3)

clear; close all;

% import code
addpath(['..' filesep 'main']);         % main
addpath(['..' filesep '45RFSOI']);      % 45rfsoi

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 1600, 470 ];
numcells    = 10;


% make object
GC = c_twoLevelGratingCell(  'discretization',   disc, ...
                            'units',            units, ...
                            'lambda',           lambda, ...
                            'domain_size',      domain, ...
                            'background_index', index_clad, ...
                            'num_cells',        numcells );

% Add a layer
height_y    = 260;
min_y       = (domain(1)-height_y)/2;
index       = 3.4;
GC          = GC.addLayer( min_y, height_y, index );

% add first rectangle
width_x     = 100;
min_x       = 0;
min_y       = min_y+20;
height_y    = 240;
index       = index_clad;
GC          = GC.addRect( min_x, min_y, width_x, height_y, index );

% add second rectangle
width_x     = 150;
min_x       = 130;
min_y       = min_y + (240-60);
height_y    = 60;
index       = index_clad;
GC          = GC.addRect( min_x, min_y, width_x, height_y, index );

% DEBUG plot the index
GC.plotIndex();

% compare bandstructures
% Pick wavelengths to sweep
period      = domain(2);
k0_max      = 10/period;
k0_min      = 0.5/period;
k0_all      = linspace( k0_min, k0_max, 100 );
lambda      = 2*pi./k0_all;

% init saving variables
all_k1      = [];   % returned k from pml type 1
all_k2      = [];   % returned k from pml type 2
all_k0      = [];
all_k_old   = [];
all_k0_old  = [];

% set simulation options
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options         = [ 1, 200, 500, 2 ];
pml_options_old     = [ 1, 200, 500, 2 ];
DEBUG               = false;

% calc bandstructure
guessk1 = 1e-5;
guessk2 = 1e-5;
guessk_old = 1e-5;
guessk  = pi/(2*period);
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
                                                               guessk, ...
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
                                                               guessk, ...
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
                                                                                   guessk, ...
                                                                                   BC, ...
                                                                                   pml_options_old );
    toc;
                 
    all_k1              = [ all_k1, k_all1.' ];           % returned k from pml type 1
    all_k2              = [ all_k2, k_all2.' ];           % returned k from pml type 2
    all_k_old           = [ all_k_old, k_all_old.' ];      % old k
    all_k0              = [ all_k0, repmat( k0_all(ii), 1, num_modes ) ];
    
%     % set new guessk's
%     guessk1         = all_k1(end);
%     guessk2         = all_k2(end);
%     guessk_old      = all_k_old(end);
    
end


% plot the bandstructure for ALL modes, new solver pml type 1 and 2
figure;
plot( real(all_k1)*period/pi, all_k0*period/pi, 'o' ); hold on;
plot( imag(all_k1)*period/pi, all_k0*period/pi, 'o' );
plot( real(all_k2)*period/pi, all_k0*period/pi, 'o' );
plot( imag(all_k2)*period/pi, all_k0*period/pi, 'o' );
plot( real(all_k_old)*period/pi, all_k0*period/pi, 'o' );
plot( imag(all_k_old)*period/pi, all_k0*period/pi, 'o' );
xlabel('ka/pi'); ylabel('k0*a/pi');
% legend('real', 'imag', 'center of bandgap');
legend('real, pm1', 'imag, pml1', 'real, pml2', 'imag, pml2', 'real, old', 'imag, old');
title('Bandstructure of all solved modes, new ver. pml 1 vs pml 2 vs. old version');
makeFigureNice();

% plot of just the old bandstructure
figure;
plot( real(all_k_old)*period/pi, all_k0*period/pi, 'o' ); hold on;
plot( imag(all_k_old)*period/pi, all_k0*period/pi, 'o' );
xlabel('ka/pi'); ylabel('k0*a/pi');
title('Bandstructure of all solved modes,old version');
makeFigureNice();

% plot of just the new bandstructure
figure;
plot( real(all_k1)*period/pi, all_k0*period/pi, 'o' ); hold on;
plot( imag(all_k1)*period/pi, all_k0*period/pi, 'o' );
plot( real(all_k2)*period/pi, all_k0*period/pi, 'o' );
plot( imag(all_k2)*period/pi, all_k0*period/pi, 'o' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real, pm1', 'imag, pml1', 'real, pml2', 'imag, pml2');
title('Bandstructure of all solved modes, new ver. pml 1 vs pml 2');
makeFigureNice();



























