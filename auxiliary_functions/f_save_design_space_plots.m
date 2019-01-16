function [] = f_save_design_space_plots( synth_obj, fig_suffix, save_plots_path )
% Saves design space plots
% Currently assumes sweep was top/bot
%
% Inputs:
%   synth_obj
%       type: c_synthTwoLeveLGrating object
%       desc: c_synthTwoLeveLGrating object, after design space has been
%             generated
%   fig_suffix
%       type: string
%       desc: string to append to figure names, if you so desire. can be empty
%             string
%   save_plots_path
%       type: string
%       desc: path to save folder

% directivity vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 10*log10(synth_obj.sweep_variables.directivities_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Directivity (dB) vs. fill factors');
figure_name = [ 'dir_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% directivity BEFORE sweeping periods vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, 10*log10(synth_obj.sweep_variables.dir_b4_period_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Directivity (dB) BEFORE PERIOD SWEEP vs. fill factors');
figure_name = [ 'dir_b4_period_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% angles vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         synth_obj.sweep_variables.angles_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Angles (deg) vs. fill factors');
figure_name = [ 'angle_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% scattering strength (imaginary) alpha vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         imag(synth_obj.sweep_variables.scatter_str_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Scattering strength (imaginary) vs. fill factors');
figure_name = [ 'scatter_str_imag_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% scattering strength alpha vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         real(synth_obj.sweep_variables.scatter_str_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Scattering strength (real) vs. fill factors');
figure_name = [ 'scatter_str_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% period vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         synth_obj.sweep_variables.periods_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title(['Period (' synth_obj.units.name ') vs. fill factors']);
figure_name = [ 'period_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% offset vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         synth_obj.sweep_variables.offsets_vs_fills );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Offset vs. fill factors');
figure_name = [ 'offsets_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

% k vs. fill
figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         real(synth_obj.sweep_variables.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Real k vs. fill factors');
figure_name = [ 'k_real_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

figure;
imagesc( synth_obj.sweep_variables.fill_tops, ...
         synth_obj.sweep_variables.fill_bots, ...
         imag(synth_obj.sweep_variables.k_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top fill factor'); ylabel('bottom fill factor');
title('Imag k vs. fill factors');
figure_name = [ 'k_imag_v_ff' fig_suffix ];
save_fig_multiformat( gcf, save_plots_path, figure_name );

end

