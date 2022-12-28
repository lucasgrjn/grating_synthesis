% authors: bohan

% Script for testing and debugging the new bloch cell class
% Geometry is that in Jelena's opt lett paper

clear; close all;

% initial settings
disc        = 10;
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

% run multiple simulations, same N, but different pmls.
% compare the k vectors for the multiple pmls
pml1 = [ 1, 200, 100, 2 ];
pml2 = [ 1, 200, 500, 2 ];
% pml3 = [ 1, 200, 1000, 2 ];

Q1 = Q.runSimulation( num_modes, BC, pml1 );
Q2 = Q.runSimulation( num_modes, BC, pml2 );
% Q3 = Q.runSimulation( num_modes, BC, pml3 );

% % sort k vectors by ascending real part
% [ ~, i_q1] = sort( real(Q1.k), 'ascend' );
% [ ~, i_q2] = sort( real(Q2.k), 'ascend' );
% [ ~, i_q3] = sort( real(Q3.k), 'ascend' );
% k_q1 = Q1.k( i_q1 );
% k_q2 = Q2.k( i_q2 );
% k_q3 = Q3.k( i_q3 );

% grab all k
k_q1 = Q1.debug.k_all;
k_q2 = Q2.debug.k_all;
% k_q3 = Q3.k;

% plot complex k and compare
k_q1_labels = {};
k_q2_labels = {};
% k_q3_labels = {};
for ii = 1:length( k_q1 )
    k_q1_labels{end+1} = [ ' ', num2str(ii) ];
end
for ii = 1:length( k_q2 )
    k_q2_labels{end+1} = [ ' ', num2str(ii) ];
end
% for ii = 1:length( k_q3 )
%     k_q3_labels{end+1} = [ ' ', num2str(ii) ];
% end

figure;
% plot q1
plot( real( k_q1 ), imag( k_q1 ), 'o' ); 
text( real( k_q1 ), imag( k_q1 ), k_q1_labels );hold on;
% plot q2
plot( real( k_q2 ), imag( k_q2 ), '+' );
text( real( k_q2 ), imag( k_q2 ), k_q2_labels );
% plot q3
% plot( real( k_q3 ), imag( k_q3 ), '*' );
% text( real( k_q3 ), imag( k_q3 ), k_q3_labels );
xlabel('real k'); ylabel('imag k');
% legend('pml type 1', 'pml type 2', 'pml type 3');
legend('pml type 1', 'pml type 2');
title('k vs. pml');
makeFigureNice();

% DEBUG plotting field
% reshape field
for ii = 1:length( k_q1 )
    % for each mode
    
    % grab field
%     field = reshape( Q1.phi_1d(:,ii), fliplr(size(Q1.N))); % dimensions (z, x) where z = direction of propagation
%     field = field.';
%     field = field./max(abs(field(:)));
    % field = fliplr( field );    % i think the field is z flipped
    field = Q1.debug.phi_all(:,:,ii);
    
    % plot field, abs
    figure;
    imagesc( Q1.x_coords, Q1.y_coords, abs( field ) );
    colorbar;
    set( gca, 'YDir', 'normal' );
    title( sprintf( 'PML TYPE 1: Field (abs) for mode %i, k real = %f', ii, real( k_q1(ii) ) ) );
    
%     % plot field, real
%     figure;
%     imagesc( Q1.x_coords, Q1.y_coords, real( field ) );
%     colorbar;
%     set( gca, 'YDir', 'normal' );
%     title( sprintf( 'PML TYPE 1: Field (real) for mode k real = %f', real( k_q1(ii) ) ) );
    
end
for ii = 1:length( k_q2 )
    % for each mode
    
    % grab field
%     field = reshape( Q2.phi_1d(:,ii), fliplr(size(Q2.N))); % dimensions (z, x) where z = direction of propagation
%     field = field.';
%     field = field./max(abs(field(:)));
    % field = fliplr( field );    % i think the field is z flipped
    field = Q2.debug.phi_all(:,:,ii);
    
    % plot field
    figure;
    imagesc( Q2.x_coords, Q2.y_coords, abs( field ) );
    colorbar;
    set( gca, 'YDir', 'normal' );
    title( sprintf( 'PML TYPE 2: Field (abs) for mode %i, k real = %f', ii, real( k_q2(ii) ) ) );
    
end
% for ii = 1:length( Q3.k )
%     % for each mode
%     
%     % grab field
%     field = reshape( Q3.phi_1d(:,ii), fliplr(size(Q3.N))); % dimensions (z, x) where z = direction of propagation
%     field = field.';
%     field = field./max(abs(field(:)));
%     % field = fliplr( field );    % i think the field is z flipped
%     
%     % plot field
%     figure;
%     imagesc( Q3.x_coords, Q3.y_coords, abs( field ) );
%     colorbar;
%     set( gca, 'YDir', 'normal' );
%     title( sprintf( 'PML TYPE 3: Field (abs) for mode k real = %f', real( Q3.k(ii) ) ) );
%     
% end















