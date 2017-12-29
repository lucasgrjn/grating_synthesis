% authors: bohan zhang
%
% script for testing reworking of FDFD complex k bloch solver
% debugging the PML
% recreating the mode from: 
%   j hu c menyuk - understanding leaky modes slab waveguide revisited
% see pg. 21 (78) of the paper for specific example used

clear; close all;

% import code
addpath(['..' filesep 'main']);                 % main
addpath(['..' filesep '45RFSOI']);              % 45rfsoi
addpath(['..' filesep 'auxiliary_functions']);  % aux functions (gui)

% initial settings
disc        = 5;
units       = 'nm';
lambda      = 1000;
k0          = 2*pi/lambda;
n0          = 1.45;
n1          = 1.39;
a           = 0.39*lambda;                              % HALF of width of core wavevguide
b           = 1.96*lambda;
outer_clad  = 1000;                                     % width of outer cladding
domain      = [ 2*outer_clad + 2*b, 10 ];
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
height_y    = 2*a;
min_y       = domain(1)/2-a;
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
num_modes   = 1;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options     = [ 1, 200, 5, 2 ];
DEBUG           = true;

% set guessk to analytical k
guessk = k0 * neff_analy;

% run new
% Phi has dimensions ny vs. nx vs. mode #
fprintf('running new solver\n');
tic;
[Phi_new, k_new, A, B] = complexk_mode_solver_2D_PML( GC.N, ...
                                                           disc, ...
                                                           k0, ...
                                                           num_modes, ...
                                                           guessk, ...
                                                           BC, ...
                                                           pml_options, ...
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

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

% display results
format shortEng                 % display in scientific notation
neff_new = k_new/k0
neff_old = k_all_old/k0

% test out the gui
f_plot_all_modes_gui( Phi_new )

% lets look at the old modes too

% generate and plot the analytical result
% constants
beta            = k0 * neff_analy;
ky              = sqrt( (k0^2)*(n0^2) - beta^2 );       % ky = kx in the paper
alpha           = sqrt( beta^2 - (k0^2)*(n1^2) );

% generate field
E_analytical    = ones( size(GC.N, 1), 1 );
x               = GC.x_coords;
y               = GC.y_coords;
y_shift         = y - y(round(end/2));                  % recenter to 0
psi             = 0;
% core section
y                                   = y_shift( abs(y_shift) < a );
E_analytical( abs(y_shift) < a )    = cos( ky * y );
% inter cladding
y                                                       = y_shift( a <= abs(y_shift) & abs(y_shift) < b );
% E_analytical( a <= abs(y_shift) & abs(y_shift) < b )    = cos( ky*a ) .* cosh( alpha*abs(y) + psi ) ./ cosh( alpha*a + psi );
E_analytical( a <= abs(y_shift) & abs(y_shift) < b )    = cos( ky*a ) .* cosh( alpha*a + psi ) ./ cosh( alpha*abs(y) + psi );
% outer cladding
y                                   = y_shift ( abs(y_shift) >= b );
% E_analytical( abs(y_shift) >= b )   = cos( ky*a ) .* cosh( alpha*b + psi ) .* exp( 1i*ky*( abs(y) - b ) ) ./ cosh( alpha*a + psi );
E_analytical( abs(y_shift) >= b )   = cos( ky*a ) .* cosh( alpha*a + psi ) .* exp( 1i*ky*( abs(y) - b ) ) ./ cosh( alpha*b + psi );

% plot analytical field xsection
figure;
plot( y_shift, real(E_analytical) );
xlabel('y, recentered'); ylabel('field, real');
title('Analytical field, real');
makeFigureNice();
% abs
figure;
plot( y_shift, abs(E_analytical) );
xlabel('y, recentered'); ylabel('field, abs');
title('Analytical field, abs');
makeFigureNice();

% plot solved field xsection
figure;
plot( y_shift, real(Phi_new(:, 1) ) );
xlabel('y, recentered'); ylabel('field, real');
title('Solved field, real');
makeFigureNice();
% plot solved field xsection
figure;
plot( y_shift, abs(Phi_new(:, 1) ) );
xlabel('y, recentered'); ylabel('field, abs');
title('Solved field, abs');
makeFigureNice();

% overplot results, normalized
figure;
plot( y_shift, abs(E_analytical)./max(abs(E_analytical)) ); hold on;
plot( y_shift, abs(Phi_new(:, 1) )./max(abs(Phi_new(:, 1) )) );
xlabel('y, recentered'); ylabel('field, abs');
legend('analytical', 'solved');
title('Analytical vs. solved field, abs');
makeFigureNice();


% -------------------------------------------------------------------------
% PML value sweeps
% -------------------------------------------------------------------------

% k vs. pml strength
pml_str = linspace( 0.1, 5, 10 );

% init saving vars
k_new_all = zeros(size(pml_str));

for ii = 1:length(pml_str)
    
    fprintf('loop %i of %i\n', ii, length(pml_str) );
    
    % set simulation options
    num_modes   = 1;
    BC          = 0;     % 0 for PEC, 1 for PMC
    % PML_options(1): PML in y direction (yes=1 or no=0)
    % PML_options(2): length of PML layer in nm
    % PML_options(3): strength of PML in the complex plane
    % PML_options(4): PML polynomial order (1, 2, 3...)
    pml_options     = [ 1, 200, pml_str(ii), 2 ];
    DEBUG           = true;

    % set guessk to analytical k
    guessk = k0 * neff_analy;
    
    % run new
    % Phi has dimensions ny vs. nx vs. mode #
    fprintf('running new solver\n');
    tic;
    [Phi_new, k_new, ~, ~] = complexk_mode_solver_2D_PML( GC.N, ...
                                                               disc, ...
                                                               k0, ...
                                                               num_modes, ...
                                                               guessk, ...
                                                               BC, ...
                                                               pml_options, ...
                                                               DEBUG );
    toc;

    % save variables
    k_new_all(ii) = k_new;
    
end

% plot effective index vs. pml strength
figure;
plot( pml_str, real(k_new_all)./k0, '-o' ); hold on;
plot( xlims, [ real(neff_analy), real(neff_analy) ], '--' );
legend('neff solver', 'analytical');
xlabel('pml strength'); ylabel('neff, real');
title('Real n_eff vs. pml_str');
makeFigureNice();
















