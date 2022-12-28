% authors: bohan

% Plot transmitted power dependence of a unit cell

clear; close all;

% dependencies
addpath(['..' filesep 'main']);
addpath(['..' filesep '45RFSOI']);

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2500, 800 ];

% directory to save data to
% unused for this script
data_dir        = '';
data_filename   = '';
data_notes      = '';

% number of parallel workers, unused
n_workers = 0;

% desired angle
optimal_angle = 15;

% coupling up/down
coupling_direction = 'down';

% make object
Q = c_synthGrating( 'discretization',   disc,       ...
                    'units',            units,      ...
                    'lambda',           lambda,     ...
                    'background_index', index_clad, ...
                    'domain_size',      domain,     ...
                    'optimal_angle',    optimal_angle,      ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers, ...
                    'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
                     );
%                     'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...


% start from waveguide mode and iterate towards the fill factor geometry i want to try
% make waveguide
waveguide = Q.h_makeGratingCell( Q.convertObjToStruct(), Q.discretization, 1.0, 1.0, 0.0 );

% run waveguide simulation
% sim settings
guess_n     = 0.7 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
guessk      = guess_n * 2*pi/Q.lambda;                                          % units rad/'units'
num_modes   = 5;
BC          = 0;                                                                % 0 = PEC
pml_options = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
% run sim
waveguide   = waveguide.runSimulation( num_modes, BC, pml_options, guessk );

% grab waveguide k
waveguide_k = waveguide.k;                                                      % units of rad/'units'     
guessk      = waveguide_k;

% DEBUG plot stuff
waveguide.plotEz_w_edges();

% calculate analytical period which would approximately phase
% match to desired output angle
k0      = Q.background_index * ( 2*pi/Q.lambda );
kx      = k0 * sin( (pi/180) * Q.optimal_angle );
period  = 2*pi/(waveguide_k- kx);                                               % units of 'units'

% snap period to discretization
guess_period    = Q.discretization * round(period/Q.discretization);


% now iterate towards the desired fill factor
des_ff  = 0.7;                                                                  % lets make the structure symmetric
ffs     = linspace( 0.95, des_ff, 10 );


% simulation settings
num_modes   = 1;
BC          = 0;                                                                % 0 = PEC
pml_options = [1, 200, 20, 2]; 
ncells      = 1;

tic;
for ii = 1:length(ffs)
    % for each fill factor
    
    fprintf('fill factor loop %i of %i...\n', ii, length(ffs) );

    % pick some periods to sweep
    sweep_periods = guess_period : disc : guess_period*1.05;
    
    % init saving variables
    angles          = zeros( size(sweep_periods) );
    k_vs_period     = zeros( size(sweep_periods) );
    GC_vs_period    = cell( size(sweep_periods) );
    
    
    for i_period = 1:length(sweep_periods)
       
        % make grating coupler object
    	GC              = Q.h_makeGratingCell( Q.convertObjToStruct(), sweep_periods(i_period), ffs(ii), ffs(ii), 0.0 );
        GC.numcells     = ncells;        % DEBUG set the number of cells to 1 for comparison
        
        % run sim
        GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

        % save angle
        if strcmp( Q.coupling_direction, 'up' )
            % coupling direction is upwards
            angles( i_period ) = GC.max_angle_up;
        else
            % coupling direction is downwards
            angles( i_period ) = GC.max_angle_down;
        end

        % update GC list
        GC_vs_period{i_period} = GC;

        % update k (units of rad/'units')
        k_vs_period(i_period)   = GC.k;
        guessk                  = GC.k;
        
    end
    
    % pick best period
    [angle_error, indx_best_period] = min( abs( Q.optimal_angle - angles ) );
    guess_period                    = sweep_periods( indx_best_period );
    guessk                          = k_vs_period( indx_best_period );
    best_GC                         = GC_vs_period{ indx_best_period };
    
    fprintf('...done\n');
    toc;
    
    
end

% save final results
best_period = guess_period;
bestk       = guessk;
    

% calc and plot power vs. x
best_GC = best_GC.power_distribution();


% plot field
best_GC.plotEz_w_edges();

% Is energy conserved? Does the radiated power + outputted power = input
% power?
total_rad_power = max(best_GC.debug.P_per_y_slice(:)) + abs(min(best_GC.debug.P_per_y_slice(:)));
P_in            = best_GC.P_in;
P_out           = best_GC.debug.P_per_x_slice(end);
P_diff          = P_in - P_out - total_rad_power;

% Assuming that the reflected power is the difference in
% input/radiated/output power, then we can calc. reflection coeff.
reflected_power_coeff = P_diff/P_in;

% what if i take the input/output powers at grid points 1/end?
P_in_grid1      = P_in * exp( 2 * best_GC.dx * imag(best_GC.k) );            % power at grid 1
P_out_gridend   = P_out * exp( -2 * best_GC.dx * imag(best_GC.k) );          % power at grid end
P_diff_1_end    = P_in_grid1 - P_out_gridend - total_rad_power;
reflected_power_coeff_1_end = P_diff_1_end/P_in_grid1;


% print results
disp('Results for not in bandgap:');
disp([ 'Total radiated power = ' num2str(total_rad_power) ]);
disp([ 'Input power, at grid point 2 = ' num2str(P_in) ]);
disp([ 'Output power, at grid point end-1 = ' num2str(P_out) ]);
disp([ 'Difference in power, grid end points - 1 = ' num2str(P_diff) ]);
disp([ 'Power reflection coeff, grid end points - 1 = ' num2str(reflected_power_coeff) ]);
disp([ 'Input power, at grid point 1 = ' num2str(P_in_grid1) ]);
disp([ 'Output power, at grid point end = ' num2str(P_out_gridend) ]);
disp([ 'Difference in power, grid end points = ' num2str(P_diff_1_end) ]);
disp([ 'Power reflection coeff, grid end points = ' num2str(reflected_power_coeff_1_end) ]);


% plot Sx and Sy distribution
figure;
imagesc( best_GC.debug.x_coords_all(2:end-1), best_GC.y_coords, (best_GC.debug.Sx) );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'S_x distribution ()' );

% plot Sx and Sy distribution
figure;
imagesc( best_GC.debug.x_coords_all, best_GC.y_coords(2:end-1), (best_GC.debug.Sy) );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'S_y distribution ()' );

% quiver/vector field plot
scale = 1.5;
figure;
quiver( best_GC.debug.x_coords_all(2:10:end-1), ...
        best_GC.y_coords(2:10:end-1), ...
        best_GC.debug.Sx(2:10:end-1, 1:10:end-1), ...
        best_GC.debug.Sy(1:10:end, 2:10:end-1), ...
        scale, 'LineWidth', 0.8, 'MaxHeadSize', 0.2 );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'Poynting vector distribution (structure outside bandgap)' );



% re-run simulation but in bandgap
% make grating coupler object
new_period                   = 650;
% best_GC_bandgap              = Q.h_makeGratingCell( Q.convertObjToStruct(), ...
%                                                     new_period, ...
%                                                     des_ff*best_GC.domain_size(2)/new_period, ...
%                                                     des_ff*best_GC.domain_size(2)/new_period, ...
%                                                     0.0 );
best_GC_bandgap              = Q.h_makeGratingCell( Q.convertObjToStruct(), ...
                                                    new_period, ...
                                                    des_ff, ...
                                                    des_ff, ...
                                                    0.0 );
best_GC_bandgap.numcells     = ncells;        % DEBUG set the number of cells to 1 for comparison
% run sim
best_GC_bandgap = best_GC_bandgap.runSimulation( num_modes, BC, pml_options, best_GC.k );

% DEBUG plot field
best_GC_bandgap.plotEz_w_edges()


% calc and plot power vs. x
best_GC_bandgap = best_GC_bandgap.power_distribution();


% calculate energy
total_rad_power = max(best_GC_bandgap.debug.P_per_y_slice(:)) + abs(min(best_GC_bandgap.debug.P_per_y_slice(:)));
P_in            = best_GC_bandgap.P_in;
P_out           = best_GC_bandgap.debug.P_per_x_slice(end);
P_diff          = P_in - P_out - total_rad_power;

% Assuming that the reflected power is the difference in
% input/radiated/output power, then we can calc. reflection coeff.
reflected_power_coeff = P_diff/P_in;

% what if i take the input/output powers at grid points 1/end?
P_in_grid1      = P_in * exp( 2 * best_GC_bandgap.dx * imag(best_GC_bandgap.k) );            % power at grid 1
P_out_gridend   = P_out * exp( -2 * best_GC_bandgap.dx * imag(best_GC_bandgap.k) );          % power at grid end
P_diff_1_end    = P_in_grid1 - P_out_gridend - total_rad_power;
reflected_power_coeff_1_end = P_diff_1_end/P_in_grid1;


% print results
disp('');
disp('Results for in bandgap:');
disp([ 'Total radiated power = ' num2str(total_rad_power) ]);
disp([ 'Input power, at grid point 2 = ' num2str(P_in) ]);
disp([ 'Output power, at grid point end-1 = ' num2str(P_out) ]);
disp([ 'Difference in power, grid end points - 1 = ' num2str(P_diff) ]);
disp([ 'Power reflection coeff, grid end points - 1 = ' num2str(reflected_power_coeff) ]);
disp([ 'Input power, at grid point 1 = ' num2str(P_in_grid1) ]);
disp([ 'Output power, at grid point end = ' num2str(P_out_gridend) ]);
disp([ 'Difference in power, grid end points = ' num2str(P_diff_1_end) ]);
disp([ 'Power reflection coeff, grid end points = ' num2str(reflected_power_coeff_1_end) ]);


% plot Sx and Sy distribution
figure;
imagesc( best_GC_bandgap.debug.x_coords_all(2:end-1), best_GC_bandgap.y_coords, (best_GC_bandgap.debug.Sx) );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'S_x distribution (structure in bandgap)' );

% plot Sx and Sy distribution
figure;
imagesc( best_GC_bandgap.debug.x_coords_all, best_GC_bandgap.y_coords(2:end-1), (best_GC_bandgap.debug.Sy) );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'S_y distribution (structure in bandgap)' );


% quiver/vector field plot
scale = 1.5;
figure;
quiver( best_GC_bandgap.debug.x_coords_all(2:10:end-1), ...
        best_GC_bandgap.y_coords(2:10:end-1), ...
        best_GC_bandgap.debug.Sx(2:10:end-1, 1:10:end-1), ...
        best_GC_bandgap.debug.Sy(1:10:end, 2:10:end-1), ...
        scale, 'LineWidth', 0.8, 'MaxHeadSize', 0.2 );
xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
set( gca, 'YDir', 'normal' );
title( 'Poynting vector distribution (structure in bandgap)' );










