% authors: bohan
% 
% script for testing the new synthesis pipeline object
%
% debugging the scattering strength calculation

clear; close all;

% add main code
addpath([ '..' filesep 'main' ]);

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1300;
index_clad  = 1.0;
domain      = [ 2000, 800 ];

% directory to save data to
% unused for this script
data_dir        = [ pwd, filesep, 'test_datasave' ];
% data_dir        = [ filesep 'project' filesep 'siphot' filesep 'bz' filesep 'gratings' filesep 'grating_synth_data' ];
data_filename   = 'test';
data_notes      = 'test sweep new dedicated function for init. grating cell';

% make the directory to save data to, if not already in existence
mkdir( data_dir );

% sweep parameters
% unused in this script
period_vec = [700, 900];
offset_vec = linspace(0, 0.3, 2);
ratio_vec  = linspace(0.7, 1.0, 1);
fill_vec   = linspace(0.5, 0.8, 1);

% number of parallel workers
n_workers = 4;

% waveguide index/thickness
waveguide_index     = [ 3.47, 3.47 ];
waveguide_thicks    = [ 100, 100 ];

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
                    'period_vec',       period_vec, ...
                    'offset_vec',       offset_vec, ...
                    'ratio_vec',        ratio_vec,  ...
                    'fill_vec',         fill_vec,   ...
                    'optimal_angle',    optimal_angle,      ...
                    'waveguide_index',  waveguide_index,    ...
                    'waveguide_thicks', waveguide_thicks,   ...
                    'coupling_direction', coupling_direction, ...
                    'data_directory',   data_dir, ...
                    'data_filename',    data_filename, ...
                    'data_notes',       data_notes, ...
                    'data_mode',        'new', ...
                    'num_par_workers',  n_workers ...
            );

% -------------------------------------------------------------------------
% loop through a range of fills and calculate the scattering strength
% -------------------------------------------------------------------------
        
% grating cell params
period  = 900;
ratio   = 1;
offset  = 0.0;

% choose fills to loop through
fills = linspace(0.5, 0.8, 4);

% init variables
power_ins   = zeros( size(fills) );
alpha_ups   = zeros( size(fills) );
alpha_dwns  = zeros( size(fills) );
p_rad_ups   = zeros( size(fills) );
p_rad_dwns  = zeros( size(fills) );
angle_ups   = zeros( size(fills) );
angle_downs = zeros( size(fills) );

for ii = 1:length(fills)
   
    fprintf( 'Loop %i of %i\n', ii, length(fills) );
    
    % make/run grating cell
    [Q, GC] = Q.testMakeGratingCell( period, fills(ii), ratio, offset );
    
    % save variables and stuff
    power_ins(ii)   = GC.P_in;
    alpha_ups(ii)   = GC.alpha_up;
    alpha_dwns(ii)  = GC.alpha_down;
    p_rad_ups(ii)   = GC.P_rad_up;
    p_rad_dwns(ii)  = GC.P_rad_down;
    angle_ups(ii)   = GC.max_angle_up;
    angle_downs(ii) = GC.max_angle_down;
    sx_up           = GC.debug.Sx_up;
    sx_down         = GC.debug.Sx_down;
    sy_up           = GC.debug.Sy_up;
    sy_down         = GC.debug.Sy_down;
    s_up            = sqrt( sx_up.^2 + sy_up(2:end-1).^2 );
    s_down          = sqrt( sx_down.^2 + sy_down(2:end-1).^2 );
    
%     % plot stuff
%     GC.plotIndex(); title( num2str(ii) );
%     GC.plotEz(); title( num2str(ii) );

    % plot poynting
    figure;
    plot( 1:length(s_up), s_up );
    xlabel('length'); ylabel('Poynting vector mag up');
    title(['poynting vector up, fill = ' num2str(fills(ii)) ]);
    makeFigureNice();
    figure;
    plot( 1:length(s_down), s_down );
    xlabel('length'); ylabel('Poynting vector mag down');
    title(['poynting vector down, fill = ' num2str(fills(ii)) ]);
    makeFigureNice();
    
end

% plot input power vs. fill
figure;
plot( fills, power_ins, '-o' );
xlabel('fill ratio'); ylabel('input power');
title('Input power vs. fill');
makeFigureNice();

% plot alphas vs. fill
figure;
plot( fills, alpha_ups, '-o', fills, alpha_dwns, '-o' );
xlabel('fill ratio'); ylabel('\alpha');
title('Alpha vs. fill');
legend('alpha up', 'alpha down');
makeFigureNice();

% plot radiated power
figure;
plot( fills, p_rad_ups, '-o', fills, p_rad_dwns, '-o' );
xlabel('fill ratio'); ylabel('power');
title('Radiated power vs. fill');
legend('up', 'down');
makeFigureNice();

% plot angle vs. fill
figure;
plot( fills, angle_ups, '-o', fills, angle_downs, '-o' );
xlabel('fill ratio'); ylabel('angle');
title('angle vs. fill');
legend('up', 'down');
makeFigureNice();

% -------------------------------------------------------------------------
% loop through a range of offsets and calculate the scattering strength
% -------------------------------------------------------------------------

% % grating cell params
% period  = 800;
% fill    = 0.8;
% ratio   = 1;
% % offset  = 0.0;
% 
% % choose offsets to loop through
% offsets = linspace(0, 0.9, 10);
% 
% % init variables
% power_ins   = zeros( size(offsets) );
% alpha_ups   = zeros( size(offsets) );
% alpha_dwns  = zeros( size(offsets) );
% p_rad_ups   = zeros( size(offsets) );
% p_rad_dwns  = zeros( size(offsets) );
% angle_ups   = zeros( size(offsets) );
% angle_downs = zeros( size(offsets) );
% 
% for ii = 1:length(offsets)
%    
%     fprintf( 'Loop %i of %i\n', ii, length(offsets) );
%     
%     % make/run grating cell
%     [Q, GC] = Q.testMakeGratingCell( period, fill, ratio, offsets(ii) );
%     
%     % save variables and stuff
%     power_ins(ii)   = GC.P_in;
%     alpha_ups(ii)   = GC.alpha_up;
%     alpha_dwns(ii)  = GC.alpha_down;
%     p_rad_ups(ii)   = GC.P_rad_up;
%     p_rad_dwns(ii)  = GC.P_rad_down;
%     angle_ups(ii)   = GC.max_angle_up;
%     angle_downs(ii) = GC.max_angle_down;
%     
%     % plot stuff
%     GC.plotIndex(); title( num2str(ii) );
%     GC.plotEz(); title( num2str(ii) );
%     
% end
% 
% % plot input power vs. offsets
% figure;
% plot( offsets, power_ins, '-o' );
% xlabel('offset'); ylabel('input power');
% title('Input power vs. offset');
% makeFigureNice();
% 
% % plot alphas vs. fill
% figure;
% plot( offsets, alpha_ups, '-o', offsets, alpha_dwns, '-o' );
% xlabel('offset'); ylabel('\alpha');
% title('Alpha vs. offset');
% legend('alpha up', 'alpha down');
% makeFigureNice();
% 
% % plot radiated power
% figure;
% plot( offsets, p_rad_ups, '-o', offsets, p_rad_dwns, '-o' );
% xlabel('offset'); ylabel('power');
% title('Radiated power vs. offset');
% legend('up', 'down');
% makeFigureNice();
% 
% % plot angle vs. fill
% figure;
% plot( offsets, angle_ups, '-o', offsets, angle_downs, '-o' );
% xlabel('offset'); ylabel('angle');
% title('angle vs. offset');
% legend('up', 'down');
% makeFigureNice();


% -------------------------------------------------------------------------
% loop through a range of periods and calculate the scattering strength
% -------------------------------------------------------------------------


% % grating cell params
% fill            = 0.8;
% ratio           = 1;
% offset          = 0.0;
% waveguide_length = 500;
% 
% % choose fills to loop through
% periods = 900:50:900;
% 
% % init variables
% power_ins   = zeros( size(periods) );
% alpha_ups   = zeros( size(periods) );
% alpha_dwns  = zeros( size(periods) );
% p_rad_ups   = zeros( size(periods) );
% p_rad_dwns  = zeros( size(periods) );
% angle_ups   = zeros( size(periods) );
% angle_downs = zeros( size(periods) );
% % sx_up       = zeros( size(periods) );
% % sx_down     = zeros( size(periods) );
% % sy_up       = zeros( size(periods) );
% % sy_down     = zeros( size(periods) );
% % s_up        = zeros( size(periods) );   % magnitude of poynting vector
% % s_down      = zeros( size(periods) );   % magnitude of poynting vector
% 
% for ii = 1:length(periods)
%    
%     fprintf( 'Loop %i of %i\n', ii, length(periods) );
%     
%     % make/run grating cell
%     [Q, GC] = Q.testMakeGratingCell( periods(ii), fill, ratio, offset );
% %     [Q, GC] = Q.testMakeGratingCell( periods(ii), waveguide_length/periods(ii), ratio, offset );    % hold wg length constant
%     
%     % save variables and stuff
%     power_ins(ii)   = GC.P_in;
%     alpha_ups(ii)   = GC.alpha_up;
%     alpha_dwns(ii)  = GC.alpha_down;
%     p_rad_ups(ii)   = GC.P_rad_up;
%     p_rad_dwns(ii)  = GC.P_rad_down;
%     angle_ups(ii)   = GC.max_angle_up;
%     angle_downs(ii) = GC.max_angle_down;
%     sx_up           = GC.debug.Sx_up;
%     sx_down         = GC.debug.Sx_down;
%     sy_up           = GC.debug.Sy_up;
%     sy_down         = GC.debug.Sy_down;
%     s_up            = sqrt( sx_up.^2 + sy_up(2:end-1).^2 );
%     s_down          = sqrt( sx_down.^2 + sy_down(2:end-1).^2 );
%     
%     % plot stuff
% %     GC.plotIndex(); title( num2str(periods(ii)) );
% %     GC.plotEz(); title( num2str(periods(ii)) );
% 
%     % plot poynting
%     figure;
%     plot( 1:length(s_up), s_up );
%     xlabel('length'); ylabel('Poynting vector mag up');
%     title(['poynting vector up, period = ' num2str(periods(ii)) ]);
%     makeFigureNice();
%     figure;
%     plot( 1:length(s_down), s_down );
%     xlabel('length'); ylabel('Poynting vector mag down');
%     title(['poynting vector down, period = ' num2str(periods(ii)) ]);
%     makeFigureNice();
%     
% end
% 
% % plot input power vs. fill
% figure;
% plot( periods, power_ins, '-o' );
% xlabel('Period (nm)'); ylabel('input power');
% title(['Input power vs. period for fill ' num2str(fill) ]);
% makeFigureNice();
% 
% % plot alphas vs. fill
% figure;
% plot( periods, alpha_ups, '-o', periods, alpha_dwns, '-o' );
% xlabel('Period (nm)'); ylabel('\alpha');
% title(['Alpha vs. period for fill ' num2str(fill) ]);
% legend('alpha up', 'alpha down');
% makeFigureNice();
% 
% % plot radiated power
% figure;
% plot( periods, p_rad_ups, '-o', periods, p_rad_dwns, '-o' );
% xlabel('Period (nm)'); ylabel('power');
% title(['Radiated power vs. period for fill ' num2str(fill) ]);
% legend('up', 'down');
% makeFigureNice();
% 
% % plot angle vs. fill
% figure;
% plot( periods, angle_ups, '-o', periods, angle_downs, '-o' );
% xlabel('Period (nm)'); ylabel('angle');
% title(['angle vs. period for fill ' num2str(fill) ]);
% legend('up', 'down');
% makeFigureNice();

























        
        
