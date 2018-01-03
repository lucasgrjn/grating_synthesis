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
outer_clad  = 920;                                      % width of outer cladding
domain      = [ 2*outer_clad + 2*b, 100 ];
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
pml_options     = [ 1, 200, 500, 2 ];
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

% ------------------
% % k vs. pml strength
% pml_str = linspace( 0.0, 50, 100 );
% 
% % init saving vars
% k_new_all = zeros(size(pml_str));
% 
% for ii = 1:length(pml_str)
%     
%     fprintf('loop %i of %i\n', ii, length(pml_str) );
%     
%     % set simulation options
%     num_modes   = 1;
%     BC          = 0;     % 0 for PEC, 1 for PMC
%     % PML_options(1): PML in y direction (yes=1 or no=0)
%     % PML_options(2): length of PML layer in nm
%     % PML_options(3): strength of PML in the complex plane
%     % PML_options(4): PML polynomial order (1, 2, 3...)
%     pml_options     = [ 1, 200, pml_str(ii), 2 ];
%     DEBUG           = false;
% 
%     % set guessk to analytical k
%     guessk = k0 * neff_analy;
%     
%     % run new
%     % Phi has dimensions ny vs. nx vs. mode #
%     fprintf('running new solver\n');
%     tic;
%     [Phi_new, k_new, ~, ~] = complexk_mode_solver_2D_PML( GC.N, ...
%                                                                disc, ...
%                                                                k0, ...
%                                                                num_modes, ...
%                                                                guessk, ...
%                                                                BC, ...
%                                                                pml_options, ...
%                                                                DEBUG );
%     toc;
% 
%     % save variables
%     k_new_all(ii) = k_new;
%     
% end
% 
% % plot effective index vs. pml strength, real
% figure;
% plot( pml_str, real(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('pml strength'); ylabel('neff, real');
% title('Real n_{eff} vs. pml strength');
% makeFigureNice();
% 
% % plot effective index vs. pml strength, imag
% figure;
% plot( pml_str, imag(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('pml strength'); ylabel('neff, imag');
% title('Imaginary n_{eff} vs. pml strength');
% makeFigureNice();
% ------------------

% % % ------------------
% % k vs. pml location
% outer_clads  = 600:20:2000;                                     % width of outer cladding
% 
% % init saving vars
% k_new_all = zeros(size(outer_clads));
% 
% for ii = 1:length(outer_clads)
%     
%     fprintf('loop %i of %i\n', ii, length(outer_clads) );
%     
%     % make new GC
%     domain = [ 2*outer_clads(ii) + 2*b, 10 ];
%     
%     % make object
%     GC = c_twoLevelGratingCell(  'discretization',   disc, ...
%                                 'units',            units, ...
%                                 'lambda',           lambda, ...
%                                 'domain_size',      domain, ...
%                                 'background_index', n1, ...
%                                 'num_cells',        numcells );
% 
%     % Add wg core
%     height_y    = 2*a;
%     min_y       = domain(1)/2-a;
%     index       = n0;
%     GC          = GC.addLayer( min_y, height_y, index );
% 
%     % Add wg bottom outer cladding
%     height_y    = outer_clad;
%     min_y       = 0;
%     index       = n0;
%     GC          = GC.addLayer( min_y, height_y, index );
% 
%     % Add wg top outer cladding
%     height_y    = outer_clad;
%     min_y       = domain(1) - outer_clad;
%     index       = n0;
%     GC          = GC.addLayer( min_y, height_y, index );
%     
%     % set simulation options
%     num_modes   = 1;
%     BC          = 0;     % 0 for PEC, 1 for PMC
%     % PML_options(1): PML in y direction (yes=1 or no=0)
%     % PML_options(2): length of PML layer in nm
%     % PML_options(3): strength of PML in the complex plane
%     % PML_options(4): PML polynomial order (1, 2, 3...)
%     pml_options     = [ 1, 200, 0.05, 2 ];
%     DEBUG           = false;
% 
%     % set guessk to analytical k
%     guessk = k0 * neff_analy;
%     
%     % run new
%     % Phi has dimensions ny vs. nx vs. mode #
%     fprintf('running new solver\n');
%     tic;
%     [Phi_new, k_new, ~, ~] = complexk_mode_solver_2D_PML( GC.N, ...
%                                                                disc, ...
%                                                                k0, ...
%                                                                num_modes, ...
%                                                                guessk, ...
%                                                                BC, ...
%                                                                pml_options, ...
%                                                                DEBUG );
%     toc;
% 
%     % save variables
%     k_new_all(ii) = k_new;
%     
% end
% 
% % plot effective index vs. pml location, real
% figure;
% plot( outer_clads, real(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('outer cladding size (nm)'); ylabel('neff, real');
% title('Real n_{eff} vs. outer cladding size/pml location');
% makeFigureNice();
% 
% % plot effective index vs. pml location, imag
% figure;
% plot( outer_clads, imag(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('outer cladding size (nm)'); ylabel('neff, imag');
% title('Imaginary n_{eff} vs. outer cladding size/pml location');
% makeFigureNice();
% % ------------------

% % ------------------
% % k vs. pml size
% pml_sizes = 50:10:500;
% 
% % init saving vars
% k_new_all = zeros(size(pml_sizes));
% 
% for ii = 1:length(pml_sizes)
%     
%     fprintf('loop %i of %i\n', ii, length(pml_sizes) );
%     
%     
%     % set simulation options
%     num_modes   = 1;
%     BC          = 0;     % 0 for PEC, 1 for PMC
%     % PML_options(1): PML in y direction (yes=1 or no=0)
%     % PML_options(2): length of PML layer in nm
%     % PML_options(3): strength of PML in the complex plane
%     % PML_options(4): PML polynomial order (1, 2, 3...)
%     pml_options     = [ 1, pml_sizes(ii), 5, 2 ];
%     DEBUG           = false;
% 
%     % set guessk to analytical k
%     guessk = k0 * neff_analy;
%     
%     % run new
%     % Phi has dimensions ny vs. nx vs. mode #
%     fprintf('running new solver\n');
%     tic;
%     [Phi_new, k_new, ~, ~] = complexk_mode_solver_2D_PML( GC.N, ...
%                                                                disc, ...
%                                                                k0, ...
%                                                                num_modes, ...
%                                                                guessk, ...
%                                                                BC, ...
%                                                                pml_options, ...
%                                                                DEBUG );
%     toc;
% 
%     % save variables
%     k_new_all(ii) = k_new;
%     
% end
% 
% % plot effective index vs. pml location, real
% figure;
% plot( pml_sizes, real(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('pml size (nm)'); ylabel('neff, real');
% title('Real n_{eff} vs. pml size');
% makeFigureNice();
% 
% % plot effective index vs. pml location, imag
% figure;
% plot( pml_sizes, imag(k_new_all)./k0, '-o' ); hold on;
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
% legend('solver', 'analytical');
% xlabel('pml size (nm)'); ylabel('neff, imag');
% title('Imaginary n_{eff} vs. opml size');
% makeFigureNice();
% % ------------------

% ------------------
% % k vs. pml strength AND ORDER
% pml_str     = linspace( 0.0, 50, 100 );
% pml_orders  = 1:6;
% 
% init saving vars
% k_new_all = zeros( length(pml_str), length(pml_orders) );       % dimensions ( pml str, pml order )
% 
% for i_order = 1:length(pml_orders)
%     for i_str = 1:length(pml_str)
% 
%         fprintf('loop %i of %i\n', i_str, length(pml_str) );
% 
%         set simulation options
%         num_modes   = 1;
%         BC          = 0;     % 0 for PEC, 1 for PMC
%         PML_options(1): PML in y direction (yes=1 or no=0)
%         PML_options(2): length of PML layer in nm
%         PML_options(3): strength of PML in the complex plane
%         PML_options(4): PML polynomial order (1, 2, 3...)
%         pml_options     = [ 1, 200, pml_str(i_str), pml_orders(i_order) ];
%         DEBUG           = false;
% 
%         set guessk to analytical k
%         guessk = k0 * neff_analy;
% 
%         run new
%         Phi has dimensions ny vs. nx vs. mode #
%         fprintf('running new solver\n');
%         tic;
%         [Phi_new, k_new, ~, ~] = complexk_mode_solver_2D_PML( GC.N, ...
%                                                                    disc, ...
%                                                                    k0, ...
%                                                                    num_modes, ...
%                                                                    guessk, ...
%                                                                    BC, ...
%                                                                    pml_options, ...
%                                                                    DEBUG );
%         toc;
% 
%         save variables
%         k_new_all(i_str, i_order) = k_new;
% 
%     end
% end
% 
% plot effective index vs. pml strength, real
% legendstrs = {};
% figure;
% for ii = 1:length( pml_orders )
%     
%     plot( pml_str, real( k_new_all(:,ii) )./k0, '-' ); hold on;
%     
%     legendstrs{end+1} = [ 'order ' num2str(ii) ];
%     
% end
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legendstrs{end+1} = 'analytical';
% legend(legendstrs);
% xlabel('pml strength'); ylabel('neff, real');
% title('Real n_{eff} vs. pml strength');
% makeFigureNice();
% 
% plot effective index vs. pml strength, imag
% figure;
% for ii = 1:length( pml_orders )
%     
%     plot( pml_str, imag( k_new_all(:,ii) )./k0, '-' ); hold on;
%     
% end
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
% legend(legendstrs);
% xlabel('pml strength'); ylabel('neff, imag');
% title('Imag n_{eff} vs. pml strength');
% makeFigureNice();
% ------------------


% ------------------
% k vs. pml location AND order
outer_clads  = 600:20:2000;                                     % width of outer cladding
% outer_clads = 4000;
pml_orders  = 1:4;

% init saving vars
k_new_all = zeros( length(outer_clads), length(pml_orders) );       % dimensions ( outer cladding, pml order )
k_old_all = k_new_all;

for i_order = 1:length(pml_orders)
    for ii = 1:length(outer_clads)

        fprintf('loop %i of %i\n', ii, length(outer_clads) );

        % make new GC
        domain = [ 2*outer_clads(ii) + 2*b, 30 ];

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
        height_y    = outer_clads(ii);
        min_y       = 0;
        index       = n0;
        GC          = GC.addLayer( min_y, height_y, index );

        % Add wg top outer cladding
        height_y    = outer_clads(ii);
        min_y       = domain(1) - outer_clads(ii);
        index       = n0;
        GC          = GC.addLayer( min_y, height_y, index );

        % set simulation options
        num_modes   = 1;
        BC          = 0;     % 0 for PEC, 1 for PMC
        % PML_options(1): PML in y direction (yes=1 or no=0)
        % PML_options(2): length of PML layer in nm
        % PML_options(3): strength of PML in the complex plane
        % PML_options(4): PML polynomial order (1, 2, 3...)
        pml_options     = [ 1, 200, 1000, pml_orders(i_order) ];
        DEBUG           = false;

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

        % save variables
        k_new_all(ii, i_order) = k_new;
        k_old_all(ii, i_order) = k_all_old;

    end
end

% solved effective index
neff_new            = k_new_all/k0;
neff_new_err_real   = 100*abs( real( neff_new - neff_analy )./real( neff_analy ) );
neff_new_err_imag   = 100*abs( imag( neff_new - neff_analy )./imag( neff_analy ) );
n_eff_old           = k_old_all/k0;
neff_old_err_real   = 100*abs( real( n_eff_old - neff_analy )./real( neff_analy ) );
neff_old_err_imag   = 100*abs( imag( n_eff_old - neff_analy )./imag( neff_analy ) );

% % plot effective index vs. pml strength, real
% legendstrs = {};
% figure;
% for ii = 1:length( pml_orders )
%     
%     plot( outer_clads, real( k_new_all(:,ii) )./k0, '-' ); hold on;
%     
%     legendstrs{end+1} = [ 'order ' num2str(ii) ];
%     
% end
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legendstrs{end+1} = 'analytical';
% legend(legendstrs);
% xlabel('outer cladding size (nm)'); ylabel('neff, real');
% title('Real n_{eff} vs. outer cladding size/pml location');
% makeFigureNice();
% 
% % plot effective index vs. pml strength, imag
% figure;
% for ii = 1:length( pml_orders )
%     
%     plot( outer_clads, imag( k_new_all(:,ii) )./k0, '-' ); hold on;
%     
% end
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
% legend(legendstrs);
% xlabel('outer cladding size (nm)'); ylabel('neff, imag');
% title('Imag n_{eff} vs. outer cladding size/pml location');
% makeFigureNice();

% plot effective index vs. pml strength, real, error
legendstrs = {};
figure;
for ii = 1:length( pml_orders )
    
    plot( outer_clads, neff_new_err_real(:,ii), '-' ); hold on;
    
    legendstrs{end+1} = [ 'order ' num2str(ii) ];
    
end
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legendstrs{end+1} = 'analytical';
legend(legendstrs);
xlabel('outer cladding size (nm)'); ylabel('% error');
title('% error b/w real n_{eff} vs. outer cladding size/pml location');
makeFigureNice();

% plot effective index vs. pml strength, imag, error
figure;
for ii = 1:length( pml_orders )
    
    plot( outer_clads, neff_new_err_imag(:,ii), '-' ); hold on;
    
end
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
legend(legendstrs);
xlabel('outer cladding size (nm)'); ylabel('% error');
title('% error b/w imag n_{eff} vs. outer cladding size/pml location');
makeFigureNice();



% plot effective index vs. pml strength, real, error, OLD
legendstrs = {};
figure;
for ii = 1:length( pml_orders )
    
    plot( outer_clads, neff_old_err_real(:,ii), '-' ); hold on;
    
    legendstrs{end+1} = [ 'order ' num2str(ii) ];
    
end
% plot( xlim, [ real(neff_analy), real(neff_analy) ], '--' );
% legendstrs{end+1} = 'analytical';
legend(legendstrs);
xlabel('outer cladding size (nm)'); ylabel('% error');
title('% error b/w real n_{eff} vs. outer cladding size/pml location, old ver.');
makeFigureNice();

% plot effective index vs. pml strength, imag, error, OLD
figure;
for ii = 1:length( pml_orders )
    
    plot( outer_clads, neff_old_err_imag(:,ii), '-' ); hold on;
    
end
% plot( xlim, [ imag(neff_analy), imag(neff_analy) ], '--' );
legend(legendstrs);
xlabel('outer cladding size (nm)'); ylabel('% error');
title('% error b/w imag n_{eff} vs. outer cladding size/pml location, old ver.');
makeFigureNice();
% ------------------






