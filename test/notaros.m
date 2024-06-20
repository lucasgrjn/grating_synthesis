% Reproduction from Bohan Zhang original code
clear; close all;

% import code
addpath(['main']);         % main
addpath(['..' filesep 'main']);         % main

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 470 ];
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

% run simulation
num_modes   = 2;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 0, 200, 500, 2 ];


% solve bandstructure
period = domain(1)
lambda_min  = 100;
lambda_max  = 5000;
k0_min      = 0.05*pi/period;
k0_max      = 3.*pi/period;
k0_all      = linspace( k0_min, k0_max, 50);

% init saving variables
all_k       = [];
all_k0      = [];
all_k_old   = [];
all_k0_old  = [];

% run solver
guessk = pi/(2*period);

x_coords    = 0:disc:period-disc;
y_coords    = 0:disc:domain(1)-disc;
% guessk = 0.000001;
for ii = 1:length(k0_all)
    
    fprintf('\nloop %i of %i\n\n', ii, length(k0_all));

    % run new
    fprintf('running new solver\n');
    tic;
    [Phi_all, k_all, A, B] = complexk_mode_solver_2D_PML( GC.N, ...
                                                           disc, ...
                                                           k0_all(ii), ...
                                                           num_modes, ...
                                                           guessk, ...
                                                           BC, ...
                                                           pml_options );
    toc;
    
    % save all k's
    all_k       = [ all_k, k_all.' ];
    all_k0      = [ all_k0, repmat( k0_all(ii), 1, length(k_all) ) ];
    
    % set new guessk
    guessk = k_all(1);
    
end

% calculate analytical bandgap properties
c           = (3e8) * 1e9;  % nm/s?
lambda_all  = 2*pi./all_k0;
k0a_pi_bg   = (2*pi/lambda)*period/pi;                            % this is where the bg should be
% gap_size_k  = k0a_pi_bg * (4/pi) * asin( abs(n2-n1)/(n2+n1) );    % size of bandgap, in k vector

% plot the bandstructure for ALL modes, new solver
figure;
plot( real(all_k)*period/pi, all_k0*period/pi, 'o' ); hold on;
plot( imag(all_k)*period/pi, all_k0*period/pi, 'o' ); hold on;
plot( xlim, [ k0a_pi_bg, k0a_pi_bg ], '--' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real', 'imag', 'center of bandgap');
title('Bandstructure of all solved modes, new ver.');
% makeFigureNice();

% plot the bandstructure for ALL modes, new solver vs wl
figure;
plot( real(all_k)*period/pi, lambda_all, 'o' ); hold on;
plot( imag(all_k)*period/pi, lambda_all, 'o' ); hold on;
xlabel('ka/pi'); ylabel('\lambda');
legend('real', 'imag');
title('Bandstructure of all solved modes, new ver.');
% makeFigureNice();

                                                                           
% reshape and sort the Phis
% when they come out raw from the modesolver, Phi_all's columns are the
% eigenvectors
% The eigenvectors are wrapped by column, then row
ny = domain(1)/disc;
nx = domain(2)/disc;

% Phi_all_half    = Phi_all( 1:end/2, : );                                        % first remove redundant bottom half
% Phi_all_reshape = reshape( Phi_all_half, ny, nx, size(Phi_all_half, 2) );       % hopefully this is dimensions y vs. x vs. mode#
% Phi_firstmode   = Phi_all_reshape( :, :, 1 );
% Phi_secondmode  = Phi_all_reshape( :, :, 2 );
Phi_firstmode   = Phi_all( :, :, 1 );
Phi_secondmode  = Phi_all( :, :, 2 );

% % do the same but with the old phi for comparison
% Phi_all_half_old    = Phi_all_old( 1:end/2, : );
% Phi_all_old_reshape = reshape( Phi_all_half_old, nx, ny, size(Phi_all_half_old, 2) );       % hopefully this is dimensions x vs. y vs. mode#
% Phi_firstmod_old    = Phi_all_old_reshape( :, :, 1 ).';                                     % gotta transpose
% 
% % x and y coords
% x_coords = 0:disc:domain(2)-disc;
% y_coords = 0:disc:domain(1)-disc;
% 
% DEBUG plot firstmode
figure;
imagesc( x_coords, y_coords, real( Phi_firstmode ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (real) for mode 1, new ver' ) );

% DEBUG plot secondmode
figure;
imagesc( x_coords, y_coords, real( Phi_secondmode ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (real) for mode 2, new ver' ) );
