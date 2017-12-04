% authors: bohan

% Script for testing and debugging the new two level grating cell class
% Calculates bandstructure
% Currently attempts to reproduce result in Jelena's opt lett paper

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1500;
index_clad  = 1.0;
domain      = [ 1600, 470 ];

% quick back of envelope calculation for period
theta   = 80 * pi/180;
neff    = 2.8;    % very approximate
kg      = 2*pi*neff/lambda;
kclad   = 2*pi*index_clad/lambda;
period_approx = (2*pi)/( kg - kclad*sin(theta) );

% run simulation
num_modes   = 20;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 500, 2 ];

% -------------------------------------------------------------------------
% Calc bandstructure
% -------------------------------------------------------------------------

% Calculating bandstructure of a grating

lambda_min  = 1000;
lambda_max  = 1e5;
n_k0        = 100;
k0_all      = linspace( 2*pi/lambda_max, 2*pi/lambda_min, n_k0 );   % equally spaced in k
lambda_all  = 2*pi./k0_all;

% make object
Q = c_twoLevelGratingCell(  'discretization', disc, ...
                            'units', units, ...
                            'lambda', lambda_all(1), ...
                            'domain_size', domain, ...
                            'background_index', index_clad );
                        
% Add a layer
height_y    = 260;
min_y       = (domain(1)-height_y)/2;
index       = 3.4;
Q           = Q.addLayer( min_y, height_y, index );

% add first rectangle
width_x     = 100;
min_x       = 0;
min_y       = min_y+20;
height_y    = 240;
index       = index_clad;
Q           = Q.addRect( min_x, min_y, width_x, height_y, index );

% add second rectangle
width_x     = 150;
min_x       = 130;
min_y       = min_y + (240-60);
height_y    = 60;
index       = index_clad;
Q           = Q.addRect( min_x, min_y, width_x, height_y, index );

% DEBUG plot the index
Q.plotIndex();

% define variables to save results in
k_all   = zeros( size(k0_all) );
phi_all = zeros( size(Q.N, 1), size(Q.N, 2), n_k0 );     % dimensions y vs x vs k0
% saves ALL the modes, even non-physical ones
k_all_all       = [];
k0_a_all_all    = [];

for ii = 1:length(lambda_all)
    % simulate once per wavelength
    
    fprintf( 'Wavelength loop %i of %i\n', ii, n_k0 );
    
    Q.lambda = lambda_all(ii);
    
    % run simulation
    Q           = Q.runSimulation( num_modes, BC, pml_options );
    k_all( ii ) = Q.k;
    phi_all(:, :, ii)   = Q.Phi;
    k_all_all           = [ k_all_all, Q.debug.k_all.' ];
    k0_a_all_all        = [ k0_a_all_all, 2*pi/lambda_all(ii)*ones(1, length( Q.debug.k_all ) ) ];
    
end

% plot bandstructure, k*a
period  = domain(2);
k0_a    = k0_all*period;
k_a     = k_all*period;
k_a_all_modes = k_all_all*period;

% Plot real bandstructure
figure;
plot( real(k_a), k0_a, 'o' );
xlabel('ka (real)'); ylabel('k_0a');
title('Bandstructure, real');
makeFigureNice();

% Plot imaginary bandstructure
figure;
plot( imag(k_a), k0_a, 'o' );
xlabel('ka (imaginary)'); ylabel('k_0a');
title('Bandstructure, imaginary');
makeFigureNice();

% plot real, ALL modes, bandstructure
figure;
plot( real(k_a_all_modes), k0_a_all_all*period, 'o' );
xlabel('ka (real)'); ylabel('k_0a');
title('Bandstructure, real, ALL MODES');
makeFigureNice();

% plot imag, ALL modes, bandstructure
figure;
plot( imag(k_a_all_modes), k0_a_all_all*period, 'o' );
xlabel('ka (imag)'); ylabel('k_0a');
title('Bandstructure, imag, ALL MODES');
makeFigureNice();

% plot field
chosen_k0 = 1;
% plot field, abs
figure;
imagesc( Q.x_coords, Q.y_coords, abs( phi_all(:,:,chosen_k0) ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Field (abs) for k0*a = %f', k0_a(chosen_k0) ) );

% Plot real vs imag ka, ALL MODES
figure;
plot( real( k_a_all_modes ), imag( k_a_all_modes ), 'o' );
xlabel('Real ka'); ylabel('imag ka');
title('real vs. imag ka, ALL MODES');
makeFigureNice();

% plot real and imag bandstructure together for pml 1
figure;
plot( real(k_a), k0_a, 'o' ); hold on;
plot( imag(k_a), k0_a, 'o' );
ylabel('k_0a');
title( 'Bandstructure, chosen modes' );
legend( 'ka real', 'ka imag' );
makeFigureNice();


% % plot bandstructure, k*a/2pi
% period  = domain(2);
% k0_a    = k0_all*period/(2*pi);
% k_a     = k_all*period/(2*pi);
% k_all_all_a2pi = k_all_all*period/(2*pi);
% 
% % Plot real bandstructure
% figure;
% plot( real(k_a), k0_a, 'o' );
% xlabel('ka/2\pi (real)'); ylabel('k_0a/2\pi');
% title('Bandstructure, real');
% makeFigureNice();
% 
% % Plot imaginary bandstructure
% figure;
% plot( imag(k_a), k0_a, 'o' );
% xlabel('ka/2\pi (imaginary)'); ylabel('k_0a/2\pi');
% title('Bandstructure, imaginary');
% makeFigureNice();
% 
% % plot real, ALL modes, bandstructure
% figure;
% plot( real(k_all_all_a2pi), k0_a_all_all, 'o' );
% xlabel('ka/2\pi (real)'); ylabel('k_0a/2\pi');
% title('Bandstructure, real, ALL MODES');
% makeFigureNice();
% 
% % plot field
% chosen_k0 = 1;
% % plot field, abs
% figure;
% imagesc( Q.x_coords, Q.y_coords, abs( phi_all(:,:,chosen_k0) ) );
% colorbar;
% set( gca, 'YDir', 'normal' );
% title( sprintf( 'Field (abs) for k0*a/2pi = %f', k0_a(chosen_k0) ) );









