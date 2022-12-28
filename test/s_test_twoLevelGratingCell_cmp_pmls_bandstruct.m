% authors: bohan

% Script for testing and debugging the new bloch cell class
% Geometry is that in Jelena's opt lett paper
%
% This version of the script compares the FULL BANDSTRUCTURE of two unit
% cells with same dielectric and different pmls

clear; close all;

% initial settings
disc        = 20;
units       = 'nm';
lambda      = 800;
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
pml_options = [ 1, 200, 100, 2 ];

% -------------------------------------------------------------------------
% Calc bandstructure
% -------------------------------------------------------------------------

% wavelength settings
lambda_min  = 1000;
lambda_max  = 1e5;
n_k0        = 100;
k0_all      = linspace( 2*pi/lambda_max, 2*pi/lambda_min, n_k0 );   % equally spaced in k
lambda_all  = 2*pi./k0_all;

% run multiple simulations, same N, but different pmls.
% compare the k vectors for the multiple pmls
pml1 = [ 1, 200, 500, 2 ];
pml2 = [ 1, 200, 1000, 2 ];

% Init a new object
% Q = c_bloch_cell( 'discretization', disc )
% Q = c_bloch_cell( 'units', units )
Q = c_twoLevelGratingCell(  'discretization', disc, ...
                            'units', units, ...
                            'lambda', lambda, ...
                            'domain_size', domain, ...
                            'background_index', index_clad )

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
k_all_1     = zeros( size(k0_all) );
phi_all_1   = zeros( size(Q.N, 1), size(Q.N, 2), n_k0 );     % dimensions y vs x vs k0
k_all_2     = zeros( size(k0_all) );
phi_all_2   = zeros( size(Q.N, 1), size(Q.N, 2), n_k0 );     % dimensions y vs x vs k0
% saves ALL the modes, even non-physical ones
k_all_all_1     = [];
k0_a_all_all_1  = [];
k_all_all_2     = [];
k0_a_all_all_2  = [];

tic;

fprintf('\nSimulating pml ver. 1:\n');
for ii = 1:length(lambda_all)
    % simulate once per wavelength
    
    fprintf( 'Wavelength loop %i of %i\n', ii, n_k0 );
    
    Q.lambda = lambda_all(ii);
    
    % run simulation
    Q                       = Q.runSimulation( num_modes, BC, pml1 );
    k_all_1( ii )           = Q.k;
    phi_all_1(:, :, ii)     = Q.Phi;
    k_all_all_1             = [ k_all_all_1, Q.debug.k_all.' ];
    k0_a_all_all_1          = [ k0_a_all_all_1, 2*pi/lambda_all(ii)*ones(1, length( Q.debug.k_all ) ) ];
    
end

fprintf('\nSimulating pml ver. 2:\n');
for ii = 1:length(lambda_all)
    % simulate once per wavelength
    
    fprintf( 'Wavelength loop %i of %i\n', ii, n_k0 );
    
    Q.lambda = lambda_all(ii);
    
    % run simulation
    Q                       = Q.runSimulation( num_modes, BC, pml2 );
    k_all_2( ii )           = Q.k;
    phi_all_2(:, :, ii)     = Q.Phi;
    k_all_all_2             = [ k_all_all_2, Q.debug.k_all.' ];
    k0_a_all_all_2          = [ k0_a_all_all_2, 2*pi/lambda_all(ii)*ones(1, length( Q.debug.k_all ) ) ];
    
end

toc;

% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

% Plot of real vs imag k
% figure;
% % plot q1
% plot( real( k_q1 ), imag( k_q1 ), 'o' ); 
% text( real( k_q1 ), imag( k_q1 ), k_q1_labels );hold on;
% % plot q2
% plot( real( k_q2 ), imag( k_q2 ), '+' );
% text( real( k_q2 ), imag( k_q2 ), k_q2_labels );
% % plot q3
% % plot( real( k_q3 ), imag( k_q3 ), '*' );
% % text( real( k_q3 ), imag( k_q3 ), k_q3_labels );
% xlabel('real k'); ylabel('imag k');
% % legend('pml type 1', 'pml type 2', 'pml type 3');
% legend('pml type 1', 'pml type 2');
% title('k vs. pml');
% makeFigureNice();

% plot bandstructure, k*a
period          = domain(2);
k0_a            = k0_all*period;
k_a_1           = k_all_1*period;
k_a_all_modes_1 = k_all_all_1*period;
k_a_2           = k_all_2*period;
k_a_all_modes_2 = k_all_all_2*period;

lambda0 = 1500;
k0_a_center = 2*pi*period/lambda0

% Plot real bandstructure
figure;
plot( real(k_a_1), k0_a, 'o' ); hold on;
plot( real(k_a_2), k0_a, 'o' );
xlabel('ka (real)'); ylabel('k_0a');
title('Bandstructure, real, chosen modes');
legend( ['pml str ' num2str(pml1(3))], ['pml str ' num2str(pml2(3))] );
makeFigureNice();

% Plot imaginary bandstructure
figure;
plot( imag(k_a_1), k0_a, 'o' ); hold on;
plot( imag(k_a_2), k0_a, 'o' );
xlabel('ka (imaginary)'); ylabel('k_0a');
title('Bandstructure, imaginary, chosen modes');
legend( ['pml str ' num2str(pml1(3))], ['pml str ' num2str(pml2(3))] );
makeFigureNice();

% plot real, ALL modes, bandstructure
figure;
plot( real(k_a_all_modes_1), k0_a_all_all_1*period, 'o' ); hold on;
plot( real(k_a_all_modes_2), k0_a_all_all_2*period, 'o' );
xlabel('ka (real)'); ylabel('k_0a');
title('Bandstructure, real, ALL MODES');
legend( ['pml str ' num2str(pml1(3))], ['pml str ' num2str(pml2(3))] );
makeFigureNice();

% plot imag, ALL modes, bandstructure
figure;
plot( imag(k_a_all_modes_1), k0_a_all_all_1*period, 'o' ); hold on;
plot( imag(k_a_all_modes_2), k0_a_all_all_2*period, 'o' );
xlabel('ka (imag)'); ylabel('k_0a');
title('Bandstructure, imag, ALL MODES');
legend( ['pml str ' num2str(pml1(3))], ['pml str ' num2str(pml2(3))] );
makeFigureNice();

% % plot field
% chosen_k0 = 1;
% % plot field, abs
% figure;
% imagesc( Q.x_coords, Q.y_coords, abs( phi_all(:,:,chosen_k0) ) );
% colorbar;
% set( gca, 'YDir', 'normal' );
% title( sprintf( 'Field (abs) for k0*a = %f', k0_a(chosen_k0) ) );
% 
% % Plot real vs imag ka, ALL MODES
% figure;
% plot( real( k_a_all_modes_1 ), imag( k_a_all_modes_1 ), 'o' );
% xlabel('Real ka'); ylabel('imag ka');
% title('real vs. imag ka, ALL MODES');
% makeFigureNice();

% let's look at the field when k0*a = 1.595, because there are 2 bands
% there. pml 1 picked a mode for band 2, and pml 2 picked a mode in band 1
chosen_k0 = 54;   % index was chosen by manual inspection
% plot field, abs, k1
figure;
imagesc( Q.x_coords, Q.y_coords, abs( phi_all_1(:,:,chosen_k0) ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Phi (abs) for pml 1, k0*a = %f', k0_a(chosen_k0) ) );
% plot field, abs, k2
figure;
imagesc( Q.x_coords, Q.y_coords, abs( phi_all_2(:,:,chosen_k0) ) );
colorbar;
set( gca, 'YDir', 'normal' );
title( sprintf( 'Phi (abs) for pml 2, k0*a = %f', k0_a(chosen_k0) ) );

% plot real and imag bandstructure together for pml 1
figure;
plot( real(k_a_1), k0_a, 'o' ); hold on;
plot( imag(k_a_1), k0_a, 'o' );
ylabel('k_0a');
title([ 'Bandstructure, chosen modes, pml str = ' num2str(pml1(3)) ]);
legend( 'ka real', 'ka imag' );
makeFigureNice();









