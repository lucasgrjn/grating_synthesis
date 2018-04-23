% authors: bohan zhang
%
% testing just the synthesis part of the new synth object
%
% let's assume a grating object has been loaded

close all;

% files to load
% lambda 1200nm angle 10 deg box 150nm
% load('D:\grating_data\2018 04 23 11 30am sorted data from scc sweep angle lambda 1200nm box 150nm mat files\angle 10 deg\synth_obj_2018_04_22_14_56_49_.mat');

% dependencies
addpath(genpath('..')); % ['..' filesep 'main']);     

% run synthesis
MFD = 10000;                                     % units nm
Q   = Q.generateFinalDesignGaussian( MFD );

% run final design in eme
Q = Q.runFinalDesignEME( MFD );

% print final results
Q.final_design

% plot final designs

% final synthesized fills
figure;
plot( 1:length(Q.bot_fill_synth), Q.bot_fill_synth, '-o' ); hold on;
plot( 1:length(Q.top_bot_fill_ratio_synth), Q.top_bot_fill_ratio_synth, '-o' );
xlabel('cell #');
legend('bottom fill ratio', 'top/bottom fill ratio');
title('Final synthesized fill factor');
makeFigureNice();
% savefig('final_fill_ratios.fig');
% saveas(gcf, 'final_fill_ratios.png');

% final synthesized fills
figure;
plot( 1:length(Q.bot_fill_synth), Q.bot_fill_synth.*Q.period_synth, '-o' ); hold on;
plot( 1:length(Q.top_bot_fill_ratio_synth), Q.top_bot_fill_ratio_synth.*Q.bot_fill_synth.*Q.period_synth, '-o' );
xlabel('cell #');
legend('bottom fill', 'top fill');
title('Final synthesized fills (nm)');
makeFigureNice();
% savefig('final_fills.fig');
% saveas(gcf, 'final_fills.png');

% final synthesized offset
figure;
plot( 1:length(Q.offset_synth), Q.offset_synth, '-o' );
xlabel('cell #');
legend('offset ratio');
title('Final synthesized offset ratio');
makeFigureNice();
% savefig('final_offsets.fig');
% saveas(gcf, 'final_offsets.png');

% final synthesized period
figure;
plot( 1:length(Q.period_synth), Q.period_synth, '-o' );
xlabel('cell #');
legend(['period, ' Q.units.name]);
title('Final synthesized period');
makeFigureNice();
% savefig('final_periods.fig');
% saveas(gcf, 'final_periods.png');

% final synthesized angles
figure;
plot( 1:length(Q.angles_synth), Q.angles_synth, '-o' );
xlabel('cell #');
legend('angle, deg');
title('Final synthesized angles');
makeFigureNice();
% savefig('final_angles.fig');
% saveas(gcf, 'final_angles.png');

% final synthesized scattering strength
figure;
plot( 1:length(Q.scatter_str_synth), Q.scatter_str_synth, '-o' );
xlabel('cell #');
legend(['scattering strength (units 1/' Q.units.name]);
title('Final synthesized scattering strength');
makeFigureNice();
% savefig('final_scattering_str.fig');
% saveas(gcf, 'final_scattering_str.png');

% % final synthesized scattering strength normalized comparison
% figure;
% plot( 1:length(Q.scatter_str_synth), Q.des_scatter_norm, '-o' ); hold on;
% plot( 1:length(Q.scatter_str_synth), Q.scatter_str_synth./abs(max(Q.scatter_str_synth(:))), '-o' );
% xlabel('cell #');
% legend('desired scatter str', 'normalized scatter str');
% title('Final synthesized scattering strength, normalized comparison');
% makeFigureNice();
% savefig('final_scattering_str_norm.fig');
% saveas(gcf, 'final_scattering_str_norm.png');

% final synthesized directivity
figure;
plot( 1:length(Q.dir_synth), 10*log10(Q.dir_synth), '-o' );
xlabel('cell #');
legend('Directivity, dB');
title('Final synthesized directivity');
makeFigureNice();
% savefig('final_dir.fig');
% saveas(gcf, 'final_dir.png');
