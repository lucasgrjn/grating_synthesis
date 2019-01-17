function [] = f_plot_save_final_design( synth_obj, save_plots_path )
% plots and saves the final design parameters vs. cell #
%
% Inputs:
%   synth_obj
%       type: c_synthTwoLevelGrating object
%       desc: c_synthTwoLevelGrating object, AFTER final design generation
%             has been run
%   save_plots_path
%       type: string
%       desc: path to save figures to

% directivity vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.dir), 10*log10(synth_obj.synthesized_design.dir), '-o' );
xlabel('unit cell #'); ylabel('directivity (dB)');
title('Directivity vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'dir_v_cell' );

% bottom fill vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.bot_fill), synth_obj.synthesized_design.bot_fill, '-o' );
xlabel('unit cell #'); ylabel('bottom fill');
title('Bottom fill vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'bot_fill_v_cell' );

% top fill vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.top_fill), synth_obj.synthesized_design.top_fill, '-o' );
xlabel('unit cell #'); ylabel('top fill');
title('Top fill vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'top_fill_v_cell' );

% top/bot fill ratio vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.top_bot_fill_ratio), synth_obj.synthesized_design.top_bot_fill_ratio, '-o' );
xlabel('unit cell #'); ylabel('top/bot fill ratio');
title('Top/bot fill ratio vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'topbot_ratio_v_cell' );

% period vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '-o' );
xlabel('unit cell #'); ylabel('period');
title('Period vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'period_v_cell' );

% offset vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.offset), synth_obj.synthesized_design.offset, '-o' );
xlabel('unit cell #'); ylabel('offset');
title('Offset vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'offset_v_cell' );

% angle vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '-o' );
xlabel('unit cell #'); ylabel('angle (deg)');
title('Angle vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'angle_v_cell' );

% scattering strength vs unit cell
figure;
plot( 1:length(synth_obj.synthesized_design.des_scatter), synth_obj.synthesized_design.des_scatter, '-o' ); hold on;
plot( 1:length(synth_obj.synthesized_design.scatter_str), synth_obj.synthesized_design.scatter_str, '-o' );
xlabel('unit cell #'); ylabel('scattering strength');
legend('Desired','Synthesized');
title('Scattering strength vs unit cell');
makeFigureNice();
save_fig_multiformat( gcf, save_plots_path, 'scatter_v_cell' );

% final index distribution
figure;
imagesc( synth_obj.synthesized_design.x_coords, synth_obj.synthesized_design.y_coords, synth_obj.synthesized_design.N );
set( gca, 'ydir', 'normal' );
colorbar;
axis image;
xlabel('nm'); ylabel('nm');
title('Final index distribution');
save_fig_multiformat( gcf, save_plots_path, 'index_distr' );

end

