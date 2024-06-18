% plotting design space for thesis

clear; close all;

% ----------------------------
% dependencies
% grating synthesis codes
addpath('C:\Users\bz\git\grating_synthesis\main');    % % on lab desktop 
addpath('C:\Users\bz\git\grating_synthesis\auxiliary_functions');    % % on lab desktop 
addpath('C:\Users\beezy\git\grating_synthesis\main'); 
addpath('C:\Users\beezy\git\grating_synthesis\auxiliary_functions');        % % on laptop     

% data
filepath = 'G:\My Drive\research\popovic group\projects\grating synthesis\data\2021 07 08 generic grating 1lvl\220108_1536_lambda1550_optangle15_dx_10_geom_fulletch';
filename = 'synth_obj.mat';

synth_obj = load( [ filepath filesep filename ] );
synth_obj = synth_obj.synth_obj;

% figure options
figopts = struct( 'font_size', 12, 'fig_size', [523,331] );
markersize = 10;

% plotting design space
figure('name', 'designspace');

yyaxis left;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, 1e3*synth_obj.sweep_variables.scatter_str_vs_fill, '.', 'markersize', markersize );
ylabel('\alpha (1/\mum)');
makeFigureNice(figopts);

yyaxis right;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, synth_obj.sweep_variables.periods_vs_fill, '.', 'markersize', markersize );
ylabel('\Lambda');
xlabel('Duty cycle');
makeFigureNice(figopts);

% pick some unit cells and plot them

% pick the cells from final design to use
dutycycles_to_plot = [ 0.25, 0.44, 0.75 ];
% chosen_cells = [ 1, 5, 8, 11 ];

% pick number of repeating unit cells to show
n_repeatunits = 4;

for i_dc = 1:length(dutycycles_to_plot)
    
    fprintf('cell %i of %i\n', i_dc, length(dutycycles_to_plot) );
    
    % find which cell to plot
    [~,i_cell] = min(abs( dutycycles_to_plot(i_dc) - synth_obj.sweep_variables.fill_ratios_to_sweep ));
    
    % grab parameters of that cell
    k           = synth_obj.sweep_variables.k_vs_fill( i_cell );
    GC          = synth_obj.sweep_variables.GC_vs_fill{ i_cell };
    GC.numcells = n_repeatunits; % set the # of repeating cells to display

    % set grating solver settings
    num_modes   = 5;
    BC          = 0;                                                % 0 = PEC
    pml_options = [1, 100, 20, 2]; 
    OPTS        = struct();
    sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, synth_obj.k0, k, OPTS );

    % plot field with waveguide overlay
    x = 1e-3 * ( 0 : GC.dx : GC.domain_size(2)*GC.numcells - GC.dx );
    y = GC.y_coords*1e-3;
    N = repmat( GC.N, 1, GC.numcells );

    figure('name', [ 'efield_dc_' strrep(num2str(dutycycles_to_plot(i_dc)), '.', 'd') ] );
    imagesc( x, y - y(end/2), real(GC.E_z) );
    colormap('redbluehilight');
    set( gca, 'ydir', 'normal' );
    axis image;
    caxis( [-1,1] .* max(abs(real(GC.E_z(:)))) );
    ylim( [-1.5, 1.5] );
    % ylim( [-1, 1] );
    % superimpose index contour
    hold all;
    contour( x, y - y(end/2), N, 1, 'color', 'k', 'LineWidth', 1.5 );
    set(gcf, 'Position', [665   391   227   302] );
    
    % plot geometry
    x = 1e-3 * ( 0 : GC.dx : GC.domain_size(2) - GC.dx );
    cmap = flipud(gray);
    figure('name', [ 'n_dc_' strrep(num2str(dutycycles_to_plot(i_dc)), '.', 'd') ], ...
            'position', [425,448,239,190] );
    imagesc( x, y - y(end/2), GC.N );
    set( gca, 'ydir', 'normal' );
    axis image;
    colormap(cmap);
    colorbar;
    caxis( [ 1, 4 ] );
    ylim( [ -1.5, 1.5 ]);
    
end

% i also want to plot just alpha and intersection with ideal alpha
w0          = 10.4;
ideal_alpha = 1/(1.4625 * w0);

figure('name', 'alpha_ideal_alpha');
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, 1e3*synth_obj.sweep_variables.scatter_str_vs_fill, '.', 'markersize', markersize ); hold on;
plot( xlim, [ ideal_alpha ideal_alpha ] );
ylabel('\alpha (1/\mum)');
xlabel('Duty cycle');
makeFigureNice(figopts);


save_all_figs('G:\My Drive\research\popovic group\thesis\figures\by section\3 2 1\working');



















