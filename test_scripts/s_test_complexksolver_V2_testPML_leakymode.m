% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% debugging the PML
% recreating the mode from: 
%   j hu c menyuk - understanding leaky modes slab waveguide revisited
% see pg. 21 (78) of the paper for specific example used

clear; close all;

% import code
addpath(['..' filesep 'main']);         % main
addpath(['..' filesep '45RFSOI']);      % 45rfsoi

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1000;
k0          = 2*pi/lambda;
n0          = 1.45;
n1          = 1.39;
a           = 0.39*lambda;
b           = 1.96*lambda;
outer_clad  = 1000;                                     % width of outer cladding
domain      = [ 2*outer_clad + a + 2*b, 200 ];
numcells    = 10;
neff_analy  = 1.4185997 + 1i*1.577e-4;


% make object
GC = c_twoLevelGratingCell(  'discretization',   disc, ...
                            'units',            units, ...
                            'lambda',           lambda, ...
                            'domain_size',      domain, ...
                            'background_index', n1, ...
                            'num_cells',        numcells );

% Add wg core
height_y    = a;
min_y       = (domain(1)-a)/2;
index       = n0;
GC          = GC.addLayer( min_y, height_y, index );

% Add wg bottom outer cladding
height_y    = outer_clad;
min_y       = 0;
index       = n0;
GC          = GC.addLayer( min_y, height_y, index );

% Add wg top outer cladding
height_y    = outer_clad;
min_y       = domain(1) - outer_clad;
index       = n0;
GC          = GC.addLayer( min_y, height_y, index );

% DEBUG plot the index
GC.plotIndex();

% -------------------------------------------------------------------------
% run sims
% -------------------------------------------------------------------------

% set simulation options
num_modes   = 10;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options     = [ 1, 200, 500, 2 ];
DEBUG           = true;

% set guessk to analytical k
guessk = k0 * neff_analy;

% run new
% Phi_all has dimensions ny vs. nx vs. mode #
fprintf('running new solver, pml type 2\n');
tic;
[Phi_newver, k_newver, A, B] = complexk_mode_solver_2D_PML( GC.N, ...
                                                           disc, ...
                                                           k0, ...
                                                           num_modes, ...
                                                           guessk, ...
                                                           BC, ...
                                                           [pml_options, 2], ...
                                                           DEBUG );
toc;

% run old
fprintf('running old solver\n');
tic;
[~, ~, Phi_all_old, k_all_old, A_old, B_old] = complexk_mode_solver_2D_PML_old( GC.N, ...
                                                                               disc, ...
                                                                               k0, ...
                                                                               num_modes, ...
                                                                               guessk, ...
                                                                               BC, ...
                                                                               pml_options );
toc;




























