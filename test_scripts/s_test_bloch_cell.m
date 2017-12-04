% authors: bohan

% Script for testing and debugging the new bloch cell class
% currently on hold

clear; close all;

% initial settings
disc    = 5e-3;
units   = 'um';
lambda  = 1.5;
indx    = 1.5;
domain  = [ 2, 5 ];

% Init a new object
% Q = c_bloch_cell( 'discretization', disc )
% Q = c_bloch_cell( 'units', units )
Q = c_bloch_cell(   'discretization', disc, ...
                    'units', units, ...
                    'lambda', lambda, ...
                    'domain_size', domain, ...
                    'background_index', indx )
                
% DEBUG plot the domain
figure;
imagesc( Q.domain.z, Q.domain.y, Q.domain.N );
xlabel( ['z (' units ')'] ); ylabel( ['y (' units ')' ]); title('DEBUG plotting dielectric');
colorbar;
set( gca, 'YDir', 'normal' );