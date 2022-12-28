% load previously saved grating synth data and plot results

clear; close all;

% add path to utility folders
addpath( ['..' filesep '..' filesep 'utility' ] );

% choose file to load
filepath = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave';
filename = '2017_11_03 19_45_32 test synth grating.mat';

% load data
load( [ filepath, filesep, filename ] );

% unpack data
v2struct( sweep_obj );
v2struct( sweep_results ); % tensors have dimensions ( fill, ratio, offset, period )

directivities(:)
angles(:)
scat_strengths(:)

% plot all the angles
figure;
plot( 1:length(angles(:)), angles(:) );
xlabel('indices'); ylabel('deg');
title('all angles');
makeFigureNice();

% let's plot all the output variables for a given period and a given fill
% ratio
period_indx = 2;
fill_indx   = 3;
directivity_mat     = squeeze( directivities(fill_indx,:,:,period_indx) );
angles_mat          = squeeze( angles(fill_indx,:,:,period_indx) );
scat_strengths_mat  = squeeze( scat_strengths(fill_indx,:,:,period_indx) );
% plot directivities
figure;
imagesc( ratio_vec, offset_vec, 10*log10(directivity_mat) );
xlabel('ratio'); ylabel('offset');
colorbar;
title( ['Directivities (dB) for period = ' num2str(period_vec(period_indx)) units ' and fill ratio = ' num2str(fill_vec(fill_indx))] );
% plot angles
figure;
imagesc( ratio_vec, offset_vec, angles_mat );
xlabel('ratio'); ylabel('offset');
colorbar;
title( ['Angles for period = ' num2str(period_vec(period_indx)) units ' and fill ratio = ' num2str(fill_vec(fill_indx))] );
% plot scat_strengths
figure;
imagesc( ratio_vec, offset_vec, 20*log10(scat_strengths_mat) );
xlabel('ratio'); ylabel('offset');
colorbar;
title( ['Scattering strength (dB) for period = ' num2str(period_vec(period_indx)) units ' and fill ratio = ' num2str(fill_vec(fill_indx))] );