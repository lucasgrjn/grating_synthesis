% authors: bohan

% Script for testing and debugging the new two level grating cell class
% sweeps offset to find best directionality

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550; %1500;
index_clad  = 1.0;
% domain      = [ 1600, 470 ];
% domain      = [ 1600, 700 ];

% quick back of envelope calculation for period
theta   = 15 * pi/180;
neff    = 2.25;    % very approximate
kg      = 2*pi*neff/lambda;
kclad   = 2*pi*index_clad/lambda;
period_approx = (2*pi)/( kg - kclad*sin(theta) );

% round the period and set the domain
period_approx_round = round( period_approx * 1e-1 )/1e-1;
period              = period_approx_round;
domain              = [ 1600, period ];

% % Init a new object
% Q = c_twoLevelGratingCell(  'discretization', disc, ...
%                             'units', units, ...
%                             'lambda', lambda, ...
%                             'domain_size', domain, ...
%                             'background_index', index_clad )


% draw two levels using two level builder function
wg_thick        = [ 80, 80 ];
wg_min_y        = [ domain(1)/2-wg_thick(1), domain(1)/2 ];
wg_indx         = [ 3.4, 3.4 ];
wgs_duty_cycles = [ 0.5, 0.5 ];

% DEBUG plot the index
% Q.plotIndex();

% run simulation
num_modes   = 20;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 500, 2 ];

% Sweep offsets and run simulations
% max_offset_top  = domain(2) - domain(2)*wgs_duty_cycles(2);
offset_bot  = domain(2) - domain(2)*wgs_duty_cycles(1);
offsets         = 0:20:domain(2);
% offsets     = 0:10:20;                                  % DEBUG 
% offsets     = [ 0, 500, 800 ];
wgs_offsets = [ offset_bot*ones(size(offsets)); offsets ];

% saving results in these variables
k_all           = [];
directivity_all = [];
p_rad_up_all    = [];
p_rad_down_all  = [];
angle_all       = [];

tic;

for ii = 1:length(offsets)
    
    fprintf('Running loop %i of %i\n', ii, length(offsets));
    
    % Init a new object
    Q = c_twoLevelGratingCell(  'discretization', disc, ...
                                'units', units, ...
                                'lambda', lambda, ...
                                'domain_size', domain, ...
                                'background_index', index_clad );
      
    % draw grating
    Q = Q.twoLevelBuilder( wg_min_y, wg_thick, wg_indx, ...
                           wgs_duty_cycles, wgs_offsets( :, ii ) );
    
    % run simulation
    Q = Q.runSimulation( num_modes, BC, pml_options );
    
    % DEBUG plot the index
%     Q.plotIndex();
    
    % save results
    k_all(end+1)            = Q.k;
    directivity_all(end+1)  = Q.directivity;
    p_rad_up_all(end+1)     = Q.P_rad_up;
    p_rad_down_all(end+1)   = Q.P_rad_down;
    angle_all(end+1)        = Q.max_angle_up;
    
    toc;
    
end

% plot directivity
figure;
plot( offsets, directivity_all, '-o' );
xlabel('Offset (um)'); ylabel('up/down directivity');
title('Up/down directivity vs. offset of top tooth from bottom tooth');
makeFigureNice();

% plot angle
figure;
plot( offsets, angle_all, '-o' );
xlabel('Offset (um)'); ylabel('output angle');
title('Output angle vs. offset of top tooth from bottom tooth');
makeFigureNice();

% plot power up
figure;
plot( offsets, p_rad_up_all, '-o' );
xlabel('Offset (um)'); ylabel('Power');
title('Power radiated up vs. offset of top tooth from bottom tooth');
makeFigureNice();

% plot power down
figure;
plot( offsets, p_rad_down_all, '-o' );
xlabel('Offset (um)'); ylabel('Power');
title('Power radiated down vs. offset of top tooth from bottom tooth');
makeFigureNice();









