function [] = f_plot_final_design( synth_obj, single_or_dual )
% plots and saves the final design parameters vs. cell #
%
% single layer curves inspired by thesis\code\2022 01 13 - plot gc results\s_plot_singlelayer_apodgc_parameters.m
% 
% Inputs:
%   synth_obj
%       type: c_synthTwoLevelGrating object
%       desc: c_synthTwoLevelGrating object, AFTER final design generation
%             has been run
%   single_or_dual
%       type: str
%       desc: 'single' or 'dual'


OPTS = struct('fig_size', [600, 600], 'grid_line_style', '-');

switch single_or_dual
    
    case 'dual'

        % Subplot version
        figure('Name', 'final_design');
        % angle vs unit cell
        subplot(3,2,1);
        plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '-o' );
        ylabel('angle (deg)');
        title('Angle vs unit cell');
        makeFigureNice( OPTS );
        % period vs unit cell
        subplot(3,2,2);
        plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '-o' );
        ylabel('period');
        title('Period vs unit cell');
        makeFigureNice( OPTS );
        % offset vs unit cell
        subplot(3,2,3);
        plot( 1:length(synth_obj.synthesized_design.offset), synth_obj.synthesized_design.offset, '-o' );
        ylabel('offset');
        title('Offset vs unit cell');
        makeFigureNice( OPTS );
        % top fill vs unit cell
        subplot(3,2,4);
        plot( 1:length(synth_obj.synthesized_design.top_fill), synth_obj.synthesized_design.top_fill, '-o' );
        ylabel('top fill');
        title('Top fill vs unit cell');
        makeFigureNice( OPTS );
        % bottom fill vs unit cell
        subplot(3,2,5);
        plot( 1:length(synth_obj.synthesized_design.bot_fill), synth_obj.synthesized_design.bot_fill, '-o' );
        ylabel('bottom fill');
        title('Bottom fill vs unit cell');
        makeFigureNice( OPTS );
        % directivity vs unit cell
        subplot(3,2,6);
        plot( 1:length(synth_obj.synthesized_design.dir), 10*log10(synth_obj.synthesized_design.dir), '-o' );
        ylabel('directivity (dB)');
        title('Directivity vs unit cell');
        makeFigureNice( OPTS );
        
    case 'single'
        
        % figure options
        OPTS = struct('fig_size', [620, 770], 'grid_line_style', '-');

        figure('Name', 'final_design_vs_cell');

        % period vs unit cell
        subplot(4,1,1);
        plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '.', 'markersize', 20 );
        xlabel('unit cell #'); ylabel('\Lambda');
        ylim( [ 0.95*min( synth_obj.synthesized_design.period ), 1.05*max(synth_obj.synthesized_design.period) ] ); 
        makeFigureNice( OPTS );

        % fill vs unit cell
        subplot(4,1,2);
        plot( 1:length(synth_obj.synthesized_design.fill), synth_obj.synthesized_design.fill, '.', 'markersize', 20 );
        xlabel('unit cell #'); ylabel('D');
        ylim( [ min(synth_obj.synthesized_design.fill) - 0.1, max(synth_obj.synthesized_design.fill) + 0.1 ] );
        makeFigureNice( OPTS );

        % angle vs unit cell
        subplot(4,1,3);
        plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '.', 'markersize', 20 );
        xlabel('unit cell #'); ylabel(['\theta (' char(176) ')']);
        ylim( [ min(synth_obj.synthesized_design.angles) - 0.5, max(synth_obj.synthesized_design.angles) + 0.5] );
        makeFigureNice( OPTS );

        % alpha vs unit cell
        subplot(4,1,4);
        plot( 1:length(synth_obj.synthesized_design.angles), 1e3*synth_obj.synthesized_design.scatter_str, '.', 'markersize', 20 );
        xlabel('unit cell #'); ylabel('\alpha (\mum^{-1})');
        ylim( [ min(1e3*synth_obj.synthesized_design.scatter_str) - max(1e3*synth_obj.synthesized_design.scatter_str)*0.2, 1.2*max(1e3*synth_obj.synthesized_design.scatter_str) ] );
        makeFigureNice( OPTS );
        
        % directivity vs unit cell
        subplot(5,1,5);
        plot( 1:length(synth_obj.synthesized_design.directivity), synth_obj.synthesized_design.directivity, '.', 'markersize', 20 );
        xlabel('unit cell #'); ylabel('Directivity (a.u)');
        ylim( [ min(synth_obj.synthesized_design.directivity) - 0.1, max(synth_obj.synthesized_design.directivity) + 0.1] );
        makeFigureNice( OPTS );

end


% % directivity vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.dir), 10*log10(synth_obj.synthesized_design.dir), '-o' );
% xlabel('unit cell #'); ylabel('directivity (dB)');
% title('Directivity vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'dir_v_cell', save_on );
% 
% % bottom fill vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.bot_fill), synth_obj.synthesized_design.bot_fill, '-o' );
% xlabel('unit cell #'); ylabel('bottom fill');
% title('Bottom fill vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'bot_fill_v_cell', save_on );
% 
% % top fill vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.top_fill), synth_obj.synthesized_design.top_fill, '-o' );
% xlabel('unit cell #'); ylabel('top fill');
% title('Top fill vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'top_fill_v_cell', save_on );
% 
% % top/bot fill ratio vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.top_bot_fill_ratio), synth_obj.synthesized_design.top_bot_fill_ratio, '-o' );
% xlabel('unit cell #'); ylabel('top/bot fill ratio');
% title('Top/bot fill ratio vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'topbot_ratio_v_cell', save_on );
% 
% % period vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '-o' );
% xlabel('unit cell #'); ylabel('period');
% title('Period vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'period_v_cell', save_on );
% 
% % offset vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.offset), synth_obj.synthesized_design.offset, '-o' );
% xlabel('unit cell #'); ylabel('offset');
% title('Offset vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'offset_v_cell', save_on );
% 
% % angle vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '-o' );
% xlabel('unit cell #'); ylabel('angle (deg)');
% title('Angle vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'angle_v_cell', save_on );

% scattering strength vs unit cell
try
    figure('Name', 'scatter_vs_cell');
    plot( 1:length(synth_obj.synthesized_design.des_scatter), synth_obj.synthesized_design.des_scatter, '-o' ); hold on;
    plot( 1:length(synth_obj.synthesized_design.scatter_str), synth_obj.synthesized_design.scatter_str, '-o' );
    xlabel('unit cell #'); ylabel('scattering strength');
    legend('Desired','Synthesized');
    title('Scattering strength vs unit cell');
    makeFigureNice();
catch
    fprintf('no scatter strength info\n');
end

% radiated power vs unit cell
try
    figure('Name', 'rad_power_ratio_vs_cell');
    plot( 1:length(synth_obj.synthesized_design.chosen_rad_power_ratio), synth_obj.synthesized_design.chosen_rad_power_ratio, '-o' );
    xlabel('unit cell #'); ylabel('radiated power ratio');
    title('radiated power ratio vs unit cell');
    makeFigureNice();
catch
    fprintf('no radiated power ratio info\n');
end
% save_fig_multiformat( gcf, save_plots_path, 'scatter_v_cell', save_on );

% final index distribution
figure('Name', 'N_synth');
imagesc( synth_obj.synthesized_design.x_coords, synth_obj.synthesized_design.y_coords, synth_obj.synthesized_design.N );
set( gca, 'ydir', 'normal' );
colorbar;
axis image;
xlabel('nm'); ylabel('nm');
title('Final index distribution');
% save_fig_multiformat( gcf, save_plots_path, 'index_distr', save_on );

% plot desired field, alpha, and realized alpha
% only if apodized design

try
    
    % figure options
    OPTS = struct('fig_size', [620, 390], 'grid_line_style', '-');

    % find where the starting point is for alpha
    % recenter xvec
    try
        xvec = 1e-3 * ( synth_obj.synthesized_design.xvec + max(synth_obj.synthesized_design.xvec) );
    catch
        xvec = 1e-3 * ( synth_obj.synthesized_design.x_coords + max(synth_obj.synthesized_design.x_coords) );
    end

    [~, indx_start] = min(abs( synth_obj.synthesized_design.alpha_des - synth_obj.synthesized_design.scatter_str(1) ) );
    x_start         = xvec(indx_start);

    alpha_des   = 1e3*synth_obj.synthesized_design.alpha_des; %1/um
    alpha_synth = 1e3*synth_obj.synthesized_design.scatter_str; %1/um

    figure('Name', 'synth_alpha_vs_x');

    plot( xvec, ...
        abs(synth_obj.synthesized_design.field_profile) .* max(alpha_des)./max(abs(synth_obj.synthesized_design.field_profile)) ); hold on;
    plot( xvec, alpha_des );
    plot( [ 0, 1e-3 * cumsum(synth_obj.synthesized_design.period(1:end-1)) ] + x_start, ...
            alpha_synth, '-o' );

    legend('Desired |E| (a.u.)', 'Desired \alpha (\mum^{-1})', 'Synth. \alpha (\mum^{-1})', ...
           'location', 'best');
    xlim([0, max(xvec)]);
    ylim([ min(alpha_des)-max(alpha_des)*0.05, max(alpha_des)*1.1]);
    xlabel('x (\mum)');
    makeFigureNice(OPTS);
    
catch
    

end

