% authors: bohan
% 
% script for testing the new synthesis pipeline object
%
% one thing that always bothered me was whether the field outputted from
% the mode solver was actually left-right flipped.

clear; close all;

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
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

% grating cell params
period  = 800;
ratio   = 1;
offset  = 0.0;
fill    = 0.80;
        
% make/run grating cell
[Q, GC] = Q.testMakeGratingCell( period, fill, ratio, offset );

% plot stuff
GC.plotIndex();
GC.plotEz();

% % trying out this edge detection business
% filt1   = 'roberts';    % roberts is the winner!
% filt2   = 'sobel';
% BW1     = edge(GC.N,filt1);   % filter
% BW2     = edge(GC.N,filt2);   % filter
% figure;
% imshowpair(BW1,BW2,'montage')
% title([filt1 ' Filter                                   ' filt2 ' Filter']);

% save variables and stuff
power_in    = GC.P_in;
alpha_up    = GC.alpha_up;
alpha_dwn   = GC.alpha_down;
p_rad_up    = GC.P_rad_up;
p_rad_dwn   = GC.P_rad_down;
angle_up    = GC.max_angle_up;
angle_down  = GC.max_angle_down;
k           = GC.k;

% use guided k to determine angle
k0 = 2*pi/lambda;
angle_from_k           = (180/pi)*asin( (real(k))/k0 )
angle_from_k_1stperiod = (180/pi)*asin( (real(k) - 2*pi/period)/k0 )

% % -------------------------------------------------------------------------
% % loop through a range of fills and calculate the scattering strength
% % -------------------------------------------------------------------------
%         
% % grating cell params
% period  = 1100;
% ratio   = 1;
% offset  = 0.0;
% 
% % choose fills to loop through
% fills = linspace(0.2, 0.9, 20);
% 
% % init variables
% power_ins   = zeros( size(fills) );
% alpha_ups   = zeros( size(fills) );
% alpha_dwns  = zeros( size(fills) );
% p_rad_ups   = zeros( size(fills) );
% p_rad_dwns  = zeros( size(fills) );
% angle_ups   = zeros( size(fills) );
% angle_downs = zeros( size(fills) );
% 
% for ii = 1:length(fills)
%    
%     fprintf( 'Loop %i of %i\n', ii, length(fills) );
%     
%     % make/run grating cell
%     [Q, GC] = Q.testMakeGratingCell( period, fills(ii), ratio, offset );
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
    



























        
        
