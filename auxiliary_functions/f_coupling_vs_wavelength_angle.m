function [ lambdas, thetas, coupling_vs_angle_lambda, overlap_vs_angle_lambda, ...
    best_coupling, best_coupling_thetas, angle_strs ] = ...
    f_coupling_vs_wavelength_angle( Ez, Hx, lambdas, n_background, center_angle, x, mfd, T )
% Calculates coupling vs. wavelength and angle
%
% Inputs:
%   Ez
%       type: double, matrix
%       desc: E field (z component) at bottom monitor
%             dimensions are x vs. wavelength
%   Hx
%       type: double, matrix
%       desc: H field (x component) at bottom monitor
%             dimensions are x vs. wavelength
%   lambda0
%       type: double, scalar
%       desc: center simulated wavelength, in meters
%   lambdas
%       type: double, array
%       desc: simulated wavelengths
%   n_background
%       type: double, scalar
%       desc: index of refraction of cladding at monitor plane
%   center_angle
%       type: double, scalar
%       desc: coupling angle which the grating was designed for
%   x
%       type: double, array
%       desc: x coordinates, in meters
%   mfd
%       type: double, scalar
%       desc: mode field diameter, in meters
%   T
%       type: double, array
%       desc: transmission from bottom field monitor vs. frequency
%
% Outputs:
%   best_coupling
%       type: double, array
%       desc: absolute best coupling at each wavelength. may not be at
%               center theta
%   best_coupling_thetas
%       type: double, array
%       desc: angles of best coupling


% if save_plots not specified as argument, assume not saving plots
if nargin < 10
    save_plots = false;
end

% squeeze data
Ez = squeeze( Ez );                                               % dimensions x vs. freq
Hx = squeeze( Hx );                                               % dimensions x vs. freq

% calculate coupling vs. angle
d0              = 0;
unit_scaling    = 1.0;
nclad           = n_background;
w0              = mfd/2; 
x_fiber         = x - x( round(end/2) );
thetas          = linspace( center_angle - 25, center_angle + 25, 101 );

% saving variables
coupling_vs_angle_lambda    = zeros( length(lambdas), length(thetas) );          % dimensions coupling vs. wavelength. vs. theta
overlap_vs_angle_lambda     = zeros( length(lambdas), length(thetas) );          % dimensions overlap vs. wavelength. vs. theta

for i_theta = 1:length(thetas)
 
    % calculate overlap vs wavelength
    field_overlap       = zeros( size(Ez) );                              % dimensions x vs. freq
    for i_wl = 1:length(lambdas)
        
        % make fiber mode
        [Ez_fiber, Hx_fiber] = f_fiberModeGaussian( w0, lambdas(i_wl), x_fiber, unit_scaling, -thetas(i_theta), d0, nclad );
        
        % calc overlap
        field_overlap(:,i_wl) = f_field_position_overlap( Ez(:,i_wl), Hx(:,i_wl), Ez_fiber, Hx_fiber );
        
    end
    
    coupling_vs_angle_lambda(:,i_theta) = T .* max( abs( field_overlap ) ).';
    overlap_vs_angle_lambda(:,i_theta)  = max( abs( field_overlap ) ).';
    
end

% plot coupling vs. angle and wavelength
% figure('Name', 'coupling_vs_ang_wl');
% imagesc( thetas, lambdas, coupling_vs_angle_lambda );
% xlabel('Angle, degrees'); ylabel('Wavelength');
% colorbar;
% set(gca, 'ydir', 'normal');
% title('Coupling vs. angle and wavelength');
% save_fig_multiformat( gcf, save_plots_path, 'coupling_vs_angle_wl', save_plots );

% plot overlap vs. angle and wavelength
% figure('Name', 'overlap_vs_ang_wl');
% imagesc( thetas, lambdas, overlap_vs_angle_lambda );
% xlabel('Angle, degrees'); ylabel('Wavelength');
% colorbar;
% set(gca, 'ydir', 'normal');
% title('Field overlap vs. angle and wavelength');
% save_fig_multiformat( gcf, save_plots_path, 'field_overlap_vs_angle_wl', save_plots );

% plot best coupling and best angle vs. wavelength
[ best_coupling, indx_best ] = max( coupling_vs_angle_lambda, [], 2 );
best_coupling_thetas         = [];
for ii = 1:size( coupling_vs_angle_lambda, 1 )
    best_coupling_thetas(ii)    = thetas(indx_best(ii));
end

 % this is for labeling the angles on the plot
angle_strs = {};
for ii = 1:length(best_coupling_thetas)
    angle_strs{ii} = [ num2str( best_coupling_thetas(ii) ) char(176) ];
end

% figure('Name', 'coup_best_ang');
% plot( lambdas, best_coupling, '-o' );
% text( lambdas, best_coupling, angle_strs, 'FontSize', 13 );
% xlabel('Wavelength'); ylabel('Coupling');
% title('Best coupling and coupling angle at each wavelength');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'best_coupling_eff', save_plots );

% % overplot best coupling and transmission
% figure('Name', 'coup_best_ang_t');
% plot( lambdas, best_coupling, '-o' ); hold on;
% text( lambdas, best_coupling, angle_strs, 'FontSize', 13 );
% plot( lambdas, T, '-o' );
% xlabel('Wavelength'); ylabel('Coupling');
% legend( 'Fiber coupling', 'Total transmission' );
% title('Best coupling and coupling angle vs. monitor transmission at each wavelength');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'best_coupling_eff_vs_T', save_plots );

% % overplot best coupling and transmission, in dB
% figure('Name', 'coup_best_ang_dB');
% plot( lambdas, 10*log10(best_coupling), '-o' ); hold on;
% text( lambdas, 10*log10(best_coupling), angle_strs, 'FontSize', 13 );
% plot( lambdas, 10*log10(T), '-o' );
% xlabel('Wavelength'); ylabel('Coupling (dB)');
% legend( 'Fiber coupling', 'Total transmission' );
% title('Best coupling and coupling angle vs. monitor transmission at each wavelength in dB');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'best_coupling_eff_vs_T_dB', save_plots );





















