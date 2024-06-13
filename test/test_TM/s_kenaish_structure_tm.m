% testing kenaish's PhC example

clear; close all;

% import code
addpath( genpath( '..' ) );

% initial settings (nm)
disc        = 10;
lambda      = 1500; %1500;
index_clad  = 1.45;
index_core  = 3.45;
y_size      = 8000;

% PhC dimensions
period      = 270;
si_width    = 1000;
gap_len     = 130;

% draw index
x_coords    = 0 : disc : period-disc;
y_coords    = 0 : disc : y_size-disc;
N           = index_clad * ones( length(y_coords), length(x_coords) );
y_high      = y_coords > y_coords(end/2) - si_width/2 - disc/2;             % waveguide
y_low       = y_coords < y_coords(end/2) + si_width/2 - disc/2;             % waveguide
% x_right     = x_coords > gap_len - disc/2;                                  % start with the gap on the left
x_left     = x_coords < x_coords(end) - gap_len - disc/2;                                  % start with the gap on the right
N( y_high & y_low, x_left ) = index_core;

% plot index
figure;
imagesc( x_coords, y_coords, N );
xlabel('x (nm)'); ylabel('y (nm)');
set( gca, 'ydir', 'normal' );
colorbar;
axis image;
title('index distribution');
          

% modesolver settings
num_modes   = 5;
BC          = 1;        % 0 for PEC, 1 for PMC
pol         = 'TM';     % 'TE' or 'TM'
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 5, 2 ];

% solve bandstructure
% lambda_min  = 100;
% lambda_max  = 5000;
lambda_all  =  1e3*[ 2.2942, 2.2333, 2.1761, 2.1222, 2.0715, 2.0236, 1.9786, 1.9362, 1.8963, 1.8589, 1.8238, 1.7910, 1.7604, 1.7321, 1.7059, 1.6821, 1.6604, 1.6411, 1.6242, 1.6096, 1.5976, 1.5863, 1.5774, 1.5720, 1.5702];
% k0_min      = 0.02*pi/period;
% k0_max      = 0.8*pi/period;
k0_min      = 2*pi/lambda_all(1);
k0_max      = 2*pi/lambda_all(end);
k0_vec      = linspace( k0_min, k0_max, length(lambda_all) );

% init saving variables
all_k       = zeros( length(k0_vec), num_modes );                               % dimensions k0 vs. mode #
all_k0      = zeros( length(k0_vec), num_modes );                               % dimensions k0 vs. mode #
% all_Phi     = zeros( length(k0_vec), ...
%                      length(y_coords), ...
%                      length(x_coords), ...
%                      num_modes );                                               % dimensions k0 vs. y vs x vs mode #
all_Phi_mode1   = zeros( length(k0_vec), length(y_coords), length(x_coords) );  % dimensions k0 vs. y vs x
all_k_mode1     = zeros( length(k0_vec), 1 ); 
all_overlaps    = zeros( length(k0_vec), num_modes );                           % dimensions k0 vs. mode #
chosen_mode_num = ones( length(k0_vec), 1 );                                    % dimensions vs. k0         
chosen_overlap  = zeros( length(k0_vec), 1 );                                   % dimensions vs. k0         

% run solver
% guessk = pi/(2*period);
% guessk = 0.0029;
% guessk = k0_min * 1.45;
% guessk = 0.25 * 2 * pi/period;
% guessk = 0.25 * 1 * pi/period;
% guessk = 0.007171;
% guessk = 0.001626 + 0.004144 * 1i;
% guessk = k0_min * 2.4;
guessk = 0;
tic;
for ii = 1:length(k0_vec)
    
    fprintf('\nloop %i of %i\n\n', ii, length(k0_vec));

    % run solver
    [Phi_solved, k_solved, ~, ~] = f_bloch_complexk_mode_solver_2D_PML( N, ...
                                                           disc, ...
                                                           k0_vec(ii), ...
                                                           num_modes, ...
                                                           guessk, ...
                                                           BC, ...
                                                           pol, ...
                                                           pml_options );
    toc;
    
%     % run solver, TE
%     [Phi_solved, k_solved, ~, ~] = f_bloch_complexk_mode_solver_2D_PML( N, ...
%                                                            disc, ...
%                                                            k0_vec(ii), ...
%                                                            num_modes, ...
%                                                            guessk, ...
%                                                            BC, ...
%                                                            'TE', ...
%                                                            pml_options );
%     toc;
    
%     k_solved( 
    
    % pick the mode with closest overlap with previous
    if ii > 1
       
        overlaps    = zeros( num_modes, 1 );
        field_1     = squeeze( all_Phi_mode1(ii-1,:,:) ) .* repmat( exp( 1i * real(all_k_mode1(ii-1)) * x_coords ), length(y_coords), 1 );
        
        for i_mode = 1:num_modes
           
            field_2             = Phi_solved(:,:,i_mode) .* repmat( exp( 1i * real(k_solved(i_mode)) * x_coords ), length(y_coords), 1 );
            overlaps(i_mode)    = f_calc_mode_overlap( field_1, field_2 );
            
        end
        
        [chosen_overlap(ii), indx_chosen_mode] = max(abs(overlaps(:)));
        
%         % idk if this matters but wrap the k into the 1st brillouin zone
%         if real( k_solved( indx_chosen_mode ) ) > pi/period
%             % k is exiting 1st brillouin zone, fold it back
%             k_solved( indx_chosen_mode ) = pi/period - mod( real(k_solved( indx_chosen_mode )), pi/period ) + 1i * imag(k_solved( indx_chosen_mode ));
%         elseif real( k_solved( indx_chosen_mode ) ) < 0
%             % k is less than 0, fold it back
%             k_solved( indx_chosen_mode ) = abs(real(k_solved( indx_chosen_mode ))) + 1i* imag(k_solved( indx_chosen_mode ));
%         end

        all_Phi_mode1( ii, :, : )   = Phi_solved(:,:,indx_chosen_mode);
        all_k_mode1( ii )           = k_solved( indx_chosen_mode );
        all_overlaps( ii, : )       = overlaps;
        chosen_mode_num(ii)         = indx_chosen_mode;
        
    else
        % first solve, pick the 1st mode
        all_Phi_mode1( ii, :, : )   = Phi_solved(:,:,1);
        all_k_mode1( ii )           = k_solved( 1 );
    end
      

    % save all k's
    all_k(ii,:)         = k_solved;
    all_k0(ii,:)        = k0_vec(ii);
    all_Phi(ii,:,:,:)   = Phi_solved;
    
    % set new guessk
    guessk = all_k_mode1( ii );
    
    % FOR DEBUGGING, plot the modes
%     f_plot_all_modes_gui(  Phi_solved, x_coords, y_coords, k_solved );
    
end

% plot the bandstructure for the chosen mode
figure;
plot( real( all_k_mode1 )*period/pi, k0_vec*period/pi, 'o' ); hold on;
plot( imag( all_k_mode1 )*period/pi, k0_vec*period/pi, 'o' ); hold on;
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ], '--' );
% % plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] + gap_size_k/2, '--' );
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] - gap_size_k/2, '--' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real', 'imag');
title('Bandstructure of chosen fundamental mode');
makeFigureNice();

% plot the bandstructure for first mode, new solver
figure;
plot( real( all_k(:,1) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
plot( imag( all_k(:,1) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ], '--' );
% % plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] + gap_size_k/2, '--' );
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] - gap_size_k/2, '--' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real', 'imag');
title('Bandstructure of first mode, new ver.');
makeFigureNice();

% plot the bandstructure for 2nd mode, new solver
figure;
plot( real( all_k(:,2) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
plot( imag( all_k(:,2) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ], '--' );
% % plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] + gap_size_k/2, '--' );
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] - gap_size_k/2, '--' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real', 'imag');
title('Bandstructure of 2nd mode, new ver.');
makeFigureNice();

% plot the bandstructure for nth mode
nth = 4;
figure;
plot( real( all_k(:,nth) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
plot( imag( all_k(:,nth) )*period/pi, k0_vec*period/pi, 'o' ); hold on;
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ], '--' );
% % plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] + gap_size_k/2, '--' );
% plot( xlim, [ k0a_pi_bg, k0a_pi_bg ] - gap_size_k/2, '--' );
xlabel('ka/pi'); ylabel('k0*a/pi');
legend('real', 'imag');
title(['Bandstructure of ' num2str(nth) 'th mode']);
makeFigureNice();


% plot gui for a chosen k0
% indx_chosen_k0 = 16;
k0_des              = 0.5455 * pi/period;
[~,indx_chosen_k0]  = min( abs( k0_vec - k0_des ) );
f_plot_all_modes_gui(  squeeze(all_Phi(indx_chosen_k0, :, :, :)), x_coords, y_coords, all_k( indx_chosen_k0) );

% plot overlap vs. iteration
figure;
plot( k0_vec, chosen_overlap, '-o' );
xlabel('k0'); ylabel('overlap');
title('Overlap vs. iteration (k0)' );
makeFigureNice();



% % run solver
% k0_cur = 0.5455;
% guessk = 0.7875*pi/period;
% num_modes = 100;
% [Phi_solved, k_solved, ~, ~] = f_bloch_complexk_mode_solver_2D_PML( N, ...
%                                                        disc, ...
%                                                        k0_cur, ...
%                                                        num_modes, ...
%                                                        guessk, ...
%                                                        BC, ...
%                                                        pol, ...
%                                                        pml_options );
% toc;










