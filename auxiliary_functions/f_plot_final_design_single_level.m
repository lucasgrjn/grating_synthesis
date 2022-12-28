function [] = f_plot_final_design_single_level( synth_obj )
% plots the final design parameters vs. cell #
% DEPRECATED
%
% Inputs:
%   synth_obj
%       type: c_synthTwoLevelGrating object
%       desc: c_synthTwoLevelGrating object, AFTER final design generation
%             has been run

OPTS = struct('fig_size', [1600, 1000]);

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
% fill vs unit cell
subplot(3,2,4);
plot( 1:length(synth_obj.synthesized_design.fill), synth_obj.synthesized_design.top_fill, '-o' );
ylabel('top fill');
title('Top fill vs unit cell');
makeFigureNice( OPTS );
% directivity vs unit cell
subplot(3,2,6);
plot( 1:length(synth_obj.synthesized_design.dir), 10*log10(synth_obj.synthesized_design.dir), '-o' );
ylabel('directivity (dB)');
title('Directivity vs unit cell');
makeFigureNice( OPTS );


% % fill vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.fill), synth_obj.synthesized_design.fill, '-o' );
% xlabel('unit cell #'); ylabel('fill');
% title('Final fill vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'fill_v_cell' );
% 
% % period vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.period), synth_obj.synthesized_design.period, '-o' );
% xlabel('unit cell #'); ylabel('period');
% title('Period vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'period_v_cell' );
% 
% % angle vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.angles), synth_obj.synthesized_design.angles, '-o' );
% xlabel('unit cell #'); ylabel('angle (deg)');
% title('Angle vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'angle_v_cell' );
% 
% % scattering strength vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.des_scatter), synth_obj.synthesized_design.des_scatter, '-o' ); hold on;
% plot( 1:length(synth_obj.synthesized_design.scatter_str), synth_obj.synthesized_design.scatter_str, '-o' );
% xlabel('unit cell #'); ylabel('scattering strength');
% legend('Desired','Synthesized');
% title('Scattering strength vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'scatter_v_cell' );
% 
% % k vs unit cell
% figure;
% plot( 1:length(synth_obj.synthesized_design.k), real( synth_obj.synthesized_design.k ), '-o' );
% xlabel('unit cell #'); ylabel('k real');
% title('k real vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'kreal_v_cell' );
% 
% figure;
% plot( 1:length(synth_obj.synthesized_design.k), imag( synth_obj.synthesized_design.k ), '-o' );
% xlabel('unit cell #'); ylabel('k imag');
% title('k real vs unit cell');
% makeFigureNice();
% save_fig_multiformat( gcf, save_plots_path, 'kimag_v_cell' );
% 
% % final index distribution
% figure;
% imagesc( synth_obj.synthesized_design.x_coords, synth_obj.synthesized_design.y_coords, synth_obj.synthesized_design.N );
% set( gca, 'ydir', 'normal' );
% colorbar;
% axis image;
% xlabel('nm'); ylabel('nm');
% title('Final index distribution');
% save_fig_multiformat( gcf, save_plots_path, 'index_distr' );




end

