function [coupling_results] = f_calc_coupling_results( synth_obj, fdtd_results )
% calculates coupling results after fdtd is ran

% angles to couple to
thetas = ( synth_obj.optimal_angle - 5 ) : 0.25 : ( synth_obj.optimal_angle + 5 );

% saving data
overlap_vs_theta_wl         = zeros( length(thetas), length(fdtd_results.lambda) );
coupling_vs_theta_wl        = zeros( size(overlap_vs_theta_wl) );
peak_overlap_vs_theta       = zeros( size(thetas) );
peak_coupling_vs_theta      = zeros( size(thetas) );
peak_overlapwl_vs_theta     = zeros( size(thetas) );
peak_couplingwl_vs_theta    = zeros( size(thetas) );
T_up_peakwl_vs_theta        = zeros( size(thetas) ); % at peak coupling wl
T_down_peakwl_vs_theta      = zeros( size(thetas) );
dir_peakwl_vs_theta         = zeros( size(thetas) );
R_peakwl_vs_theta           = zeros( size(thetas) );
coupling_1db_bw_vs_theta    = zeros( size(thetas) );
overlap_1db_bw_vs_theta     = zeros( size(thetas) );

directionality_vs_wl = fdtd_results.T_up./fdtd_results.T_down;
if strcmp( synth_obj.coupling_direction, 'down' )
    directionality_vs_wl = 1./directionality_vs_wl;
end

% calculate coupling results for each angle
for i_theta = 1:length(thetas)

    switch synth_obj.coupling_direction

        % depending on coupling direction grab the right data
        case 'up'
            Ez      = fdtd_results.Ez_up;
            Hx      = fdtd_results.Hx_up;
            theta   = -thetas(i_theta);
            T       = fdtd_results.T_up.';

        case 'down'
            Ez      = fdtd_results.Ez_down;
            Hx      = -fdtd_results.Hx_down;
            theta   = -thetas(i_theta);
            T       = fdtd_results.T_down.';

    end % end coupling direction switch
    
    % calculate overlap
    [overlap_vs_theta_wl(i_theta,:)] = f_overlap_1d( Ez, Hx, 1e-9 * synth_obj.synthesized_design.MFD/2, ...
                          fdtd_results.lambda, theta, ...
                          fdtd_results.x, synth_obj.background_index );

    % calculate coupling
    coupling_vs_theta_wl(i_theta,:) = overlap_vs_theta_wl(i_theta,:) .* T;

    % find value, wavelength, and 1dB bw of peak coupling
    [ peak_coupling_vs_theta(i_theta),...
      peak_couplingwl_vs_theta(i_theta),...
      indx_peak, ...
      coupling_1db_bw_vs_theta(i_theta)] = f_getpeak_calc_1db_bw( fdtd_results.lambda, coupling_vs_theta_wl(i_theta,:) );

    % get other values that we care about at peak coupling wavelength
    T_up_peakwl_vs_theta(i_theta)   = fdtd_results.T_up(indx_peak);
    T_down_peakwl_vs_theta(i_theta) = fdtd_results.T_down(indx_peak);
    dir_peakwl_vs_theta(i_theta)    = fdtd_results.T_up(indx_peak)./fdtd_results.T_down(indx_peak);
    if strcmp( synth_obj.coupling_direction, 'down' )
        dir_peakwl_vs_theta(i_theta) = 1./dir_peakwl_vs_theta(i_theta);
    end
    R_peakwl_vs_theta(i_theta)      = fdtd_results.R(indx_peak);

    
    % find value, wavelength, and 1dB bw of peak overlap
    [ peak_overlap_vs_theta(i_theta),...
      peak_overlapwl_vs_theta(i_theta),...
      ~, ...
      overlap_1db_bw_vs_theta(i_theta)] = f_getpeak_calc_1db_bw( fdtd_results.lambda, overlap_vs_theta_wl(i_theta,:) );

end % end theta for loop

% grab some specific results that i care about

% results at best coupling angle for the center wavelength
[~, indx_center_wl] = min(abs(fdtd_results.lambda - 1e-9*synth_obj.lambda));
center_wl           = fdtd_results.lambda(indx_center_wl);

% T and R at the center wavelength
T_up_centerwl    = fdtd_results.T_up(indx_center_wl);
T_down_centerwl  = fdtd_results.T_down(indx_center_wl);
dir_centerwl     = fdtd_results.T_up(indx_center_wl)./fdtd_results.T_down(indx_center_wl);
if strcmp( synth_obj.coupling_direction, 'down' )
    dir_centerwl = 1./dir_centerwl;
end
R_centerwl = fdtd_results.R(indx_center_wl);

% overlap, coupling, and bandwidth at best coupling angle for the center wl
coupling_vs_theta_centerwl       = coupling_vs_theta_wl(:,indx_center_wl);
[best_coupling_centerwl, indx_ang_bestcoup_centerwl] = max(coupling_vs_theta_centerwl); % get best coupling and index of angle at center wavelength
overlap_atbestcoupling_centerwl  = overlap_vs_theta_wl( indx_ang_bestcoup_centerwl, indx_center_wl );
best_coupling_ang_centerwl       = thetas(indx_ang_bestcoup_centerwl);
overlap_1db_bw_bestang_centerwl  = overlap_1db_bw_vs_theta(indx_ang_bestcoup_centerwl);
coupling_1db_bw_bestang_centerwl = coupling_1db_bw_vs_theta(indx_ang_bestcoup_centerwl);

% results at the best coupling angle and wavelength
best_coupling_alldata         = max(coupling_vs_theta_wl(:)); % best coupling across all angles and wavelength
[indx_besttheta, indx_bestwl] = find( coupling_vs_theta_wl == best_coupling_alldata );
bestwl_alldata                = fdtd_results.lambda(indx_bestwl);
bestang_alldata               = thetas(indx_besttheta);
overlap_atbestcoup_alldata    = overlap_vs_theta_wl( indx_besttheta, indx_bestwl );
overlap_1db_bw_alldata        = overlap_1db_bw_vs_theta(indx_besttheta);
coupling_1db_bw_alldata       = coupling_1db_bw_vs_theta(indx_besttheta);
T_up_bestwl_alldata           = fdtd_results.T_up(indx_bestwl);
T_down_bestwl_alldata         = fdtd_results.T_down(indx_bestwl);
dir_bestwl_alldata            = fdtd_results.T_up(indx_bestwl)./fdtd_results.T_down(indx_bestwl);
if strcmp( synth_obj.coupling_direction, 'down' )
    dir_bestwl_alldata = 1./dir_bestwl_alldata;
end
R_bestwl_alldata = fdtd_results.R(indx_bestwl);

% results at the best coupling wavelength for the intended design angle
[ ~, indx_designang ]           = min(abs(synth_obj.optimal_angle-thetas));
design_ang                      = synth_obj.optimal_angle;
[coupling_designang, indx_wl_bestcoup_designang] = max(coupling_vs_theta_wl(indx_designang,:));
wl_bestcoup_designang           = fdtd_results.lambda(indx_wl_bestcoup_designang);
overlap_atbestcoup_designang    = overlap_vs_theta_wl(indx_designang, indx_wl_bestcoup_designang);
overlap_1db_bw_designang        = overlap_1db_bw_vs_theta(indx_designang);
coupling_1db_bw_designang       = coupling_1db_bw_vs_theta(indx_designang);
T_up_bestcoup_designang         = fdtd_results.T_up(indx_wl_bestcoup_designang);
T_down_bestcoup_designang       = fdtd_results.T_down(indx_wl_bestcoup_designang);
dir_bestcoup_designang          = fdtd_results.T_up(indx_wl_bestcoup_designang)./fdtd_results.T_down(indx_wl_bestcoup_designang);
if strcmp( synth_obj.coupling_direction, 'down' )
    dir_bestcoup_designang = 1./dir_bestcoup_designang;
end
R_bestcoup_designang = fdtd_results.R(indx_wl_bestcoup_designang);

% find 1dB bw of T
[ peak_T,...
  peak_T_wl,...
  indx_peak, ...
  T_1db_bw ] = f_getpeak_calc_1db_bw( fdtd_results.lambda, T );

% find where reflection increases to 5%
% [~, ~, ~, R_5perc_bw] = f_getpeak_calc_bw( wl, val_vs_wl, thresh )

% results to return
coupling_results = struct( ...
    'thetas', thetas, ...
    'directionality_vs_wl', directionality_vs_wl, ...
    'overlap_vs_theta_wl', overlap_vs_theta_wl, ...
    'coupling_vs_theta_wl', coupling_vs_theta_wl, ...
    'peak_overlap_vs_theta', peak_overlap_vs_theta, ...
    'peak_coupling_vs_theta', peak_coupling_vs_theta, ...
    'peak_overlapwl_vs_theta', 1e9*peak_overlapwl_vs_theta, ...
    'peak_couplingwl_vs_theta', 1e9*peak_couplingwl_vs_theta, ...
    'peak_T', peak_T, ...
    'peak_T_wl', 1e9*peak_T_wl, ...
    'T_1db_bw', 1e9*T_1db_bw, ...
    'T_up_peakwl_vs_theta', T_up_peakwl_vs_theta, ...
    'T_down_peakwl_vs_theta', T_down_peakwl_vs_theta, ...
    'dir_peakwl_vs_theta', dir_peakwl_vs_theta, ...
    'R_peakwl_vs_theta', R_peakwl_vs_theta, ...
    'coupling_1db_bw_vs_theta', 1e9*coupling_1db_bw_vs_theta, ...
    'overlap_1db_bw_vs_theta', 1e9*overlap_1db_bw_vs_theta, ...
    'center_wl', 1e9*center_wl, ... 
    'best_coupling_ang_centerwl', best_coupling_ang_centerwl, ...
    'T_up_centerwl', T_up_centerwl, ...
    'T_up_centerwl_dB', 10*log10(T_up_centerwl), ...
    'T_down_centerwl', T_down_centerwl, ...
    'T_down_centerwl_dB', 10*log10(T_down_centerwl), ...
    'dir_centerwl', dir_centerwl, ...
    'dir_centerwl_dB', 10*log10(dir_centerwl), ...
    'R_centerwl', R_centerwl, ...
    'R_centerwl_IL_dB', 10*log10(1 - abs(R_centerwl)), ...
    'best_coupling_centerwl', best_coupling_centerwl, ...
    'best_coupling_centerwl_dB', 10*log10(best_coupling_centerwl), ...
    'overlap_atbestcoupling_centerwl', overlap_atbestcoupling_centerwl, ...
    'overlap_atbestcoupling_centerwl_dB', 10*log10(overlap_atbestcoupling_centerwl), ...
    'overlap_1db_bw_bestang_centerwl', 1e9*overlap_1db_bw_bestang_centerwl, ...
    'coupling_1db_bw_bestang_centerwl', 1e9*coupling_1db_bw_bestang_centerwl, ...
    'bestwl_alldata', 1e9*bestwl_alldata, ...
    'bestang_alldata', bestang_alldata, ...
    'T_up_bestwl_alldata', T_up_bestwl_alldata, ...
    'T_up_bestwl_alldata_dB', 10*log10(T_up_bestwl_alldata), ...
    'T_down_bestwl_alldata', T_down_bestwl_alldata, ...
    'T_down_bestwl_alldata_dB', 10*log10(T_down_bestwl_alldata), ...
    'dir_bestwl_alldata', dir_bestwl_alldata, ...
    'dir_bestwl_alldata_dB', 10*log10(dir_bestwl_alldata), ...
    'R_bestwl_alldata', R_bestwl_alldata, ...
    'R_bestwl_alldata_IL_dB', 10*log10(1-abs(R_bestwl_alldata)), ...
    'best_coupling_alldata', best_coupling_alldata, ...
    'best_coupling_alldata_dB', 10*log10(best_coupling_alldata), ...
    'overlap_atbestcoup_alldata', overlap_atbestcoup_alldata, ...
    'overlap_atbestcoup_alldata_dB', 10*log10(overlap_atbestcoup_alldata), ...
    'overlap_1db_bw_alldata', 1e9*overlap_1db_bw_alldata, ...
    'coupling_1db_bw_alldata', 1e9*coupling_1db_bw_alldata, ...
    'design_ang', design_ang, ...
    'wl_bestcoup_designang', 1e9*wl_bestcoup_designang, ...
    'T_up_bestcoup_designang', T_up_bestcoup_designang, ...
    'T_up_bestcoup_designang_dB', 10*log10(T_up_bestcoup_designang), ...
    'T_down_bestcoup_designang', T_down_bestcoup_designang, ...
    'T_down_bestcoup_designang_dB', 10*log10(T_down_bestcoup_designang), ...
    'dir_bestcoup_designang', dir_bestcoup_designang, ...
    'dir_bestcoup_designang_dB', 10*log10(dir_bestcoup_designang), ...
    'R_bestcoup_designang', R_bestcoup_designang, ...
    'R_bestcoup_designang_IL_dB', 10*log10(1-abs(R_bestcoup_designang)), ...
    'coupling_designang', coupling_designang, ...
    'coupling_designang_dB', 10*log10(coupling_designang), ...
    'overlap_atbestcoup_designang', overlap_atbestcoup_designang, ...
    'overlap_atbestcoup_designang_dB', 10*log10(overlap_atbestcoup_designang), ...
    'overlap_1db_bw_designang', 1e9*overlap_1db_bw_designang, ...
    'coupling_1db_bw_designang', 1e9*coupling_1db_bw_designang, ...
    'indx_designang', indx_designang, ...
    'indx_ang_bestcoup_centerwl', indx_ang_bestcoup_centerwl, ...
    'indx_besttheta', indx_besttheta, ... 
    'indx_bestwl', indx_bestwl ); 

end




% helper function for getting bw
function [peak_val, peak_wl, indx_peak, bw_1db] = f_getpeak_calc_1db_bw( wl, val_vs_wl )
% first finds the value, index, and wavelength of maximum
% then calculates 1db bw
% 
% wl = wavelength
% val_vs_wl = function in amplitude not db

[peak_val, indx_peak]   = max(val_vs_wl);
peak_wl                 = wl(indx_peak);

% get wavelength bounds
try
    wl_bound1 = interp1( 10*log10( val_vs_wl(1:indx_peak) ), ...
                               wl(1:indx_peak), ...
                               10*log10( peak_val ) - 1 );
catch
    % i think error gets thrown when peak is at edge of sim domain
    wl_bound1 = wl(1);
end
try
    wl_bound2 = interp1( 10*log10( val_vs_wl(indx_peak:end) ), ...
                               wl(indx_peak:end), ...
                               10*log10( peak_val ) - 1 );
catch
    % i think error gets thrown when peak is at edge of sim domain
    wl_bound2 = wl(end);
end
                       
% check if wavelengths are out of bands and if so, snap to edge of
% simulated bandwidth
if isnan(wl_bound1)
    wl_bound1 = wl(1);
end
if isnan(wl_bound2)
    wl_bound2 = wl(end);
end
                       
bw_1db = abs(wl_bound2 - wl_bound1);

end

% helper function for getting bw
function [peak_val, peak_wl, indx_peak, bw_1db] = f_getpeak_calc_bw( wl, val_vs_wl, thresh )
% first finds the value, index, and wavelength of maximum
% then calculates 1db bw
% 
% wl = wavelength
% val_vs_wl = function in amplitude not db
% thresh = threshold (down from max) to calc bw for, in dB

[peak_val, indx_peak]   = max(val_vs_wl);
peak_wl                 = wl(indx_peak);

% get wavelength bounds
try
    wl_bound1 = interp1( 10*log10( val_vs_wl(1:indx_peak) ), ...
                               wl(1:indx_peak), ...
                               10*log10( peak_val ) - thresh );
catch
    % i think error gets thrown when peak is at edge of sim domain
    wl_bound1 = wl(1);
end
try
    wl_bound2 = interp1( 10*log10( val_vs_wl(indx_peak:end) ), ...
                               wl(indx_peak:end), ...
                               10*log10( peak_val ) - thresh );
catch
    % i think error gets thrown when peak is at edge of sim domain
    wl_bound2 = wl(end);
end
                       
% check if wavelengths are out of bands and if so, snap to edge of
% simulated bandwidth
if isnan(wl_bound1)
    wl_bound1 = wl(1);
end
if isnan(wl_bound2)
    wl_bound2 = wl(end);
end
                       
bw_1db = abs(wl_bound2 - wl_bound1);

end




















