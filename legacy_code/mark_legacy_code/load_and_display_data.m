% loading and displaying saved solver data

clear; clc; close all;

% ------------------------------------------------------------------------|
% Loading
% ------------------------------------------------------------------------|

% pick datafile to load
PathName = 'C:\Users\bozh7002\Google Drive\research\popovic group\data\gratings\';  % plab desktop
FileName = '2017_01_23 16_38_35 data best mode.mat';

% load data
load( [PathName,FileName] );

n_runs = length(all_data.runs);

% ------------------------------------------------------------------------|
% Plotting
% ------------------------------------------------------------------------|

h_env = figure;
h_tot = figure;



% plot all fields
for i_run = 1:n_runs

    % get current run data
    cur_data    = all_data.runs(i_run);
    field_all   = cur_data.field;         % dimensions are (x, z, mode#)
    x_vec       = cur_data.x_vec;
    z_vec       = cur_data.z_vec;
    k           = cur_data.k;
    k0          = cur_data.k0;
    guessk      = cur_data.guessk;
    dx          = cur_data.dx;
    lambda0     = cur_data.lambda0;
    try
        b = 2*pi/cur_data.bloch_period;
    catch
        b = 2*pi/650;       % no bloch period found in data, default to 650nm
        fprintf('No bloch_period found, defaulting to 650nm\n');
    end
    
    % plot field, every mode
    for i_mode = 1:length(k)
        
        % normalize current field
        cur_field = field_all(:,:,i_mode);
        cur_field = cur_field./max(abs(cur_field(:)));
        
        % add phase to current field
        k_phase             = exp(1i*k(i_mode).*z_vec);
        cur_field_wphase    = cur_field.*repmat(k_phase, length(x_vec), 1);
        
        % extend field across several periods, include phase
        n_periods       = 20;
        cur_field_rep   = repmat(cur_field,1,n_periods);
        z_vec_long      = linspace(0, n_periods*max(z_vec), size(cur_field_rep,2) );
        k_phase         = exp(1i*k(i_mode).*z_vec_long);
        cur_field_rep   = cur_field_rep.*repmat(k_phase,length(x_vec),1);
        
        % plot field, no phase
        figure(h_env); set(gcf,'Position',[50 50 1700 900]);
        subplot(1, 3, 1);
        imagesc(z_vec, x_vec, imag(cur_field) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Bloch envelope (imag), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        subplot(1, 3, 2);
        imagesc(z_vec, x_vec, real(cur_field) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Bloch envelope (real), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        subplot(1, 3, 3);
        imagesc(z_vec, x_vec, abs(cur_field) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Bloch envelope (amplitude), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        % plot field, with phase, across several periods
        figure(h_tot); set(gcf,'Position',[50 50 1700 900]);
        subplot(1, 3, 1);
        imagesc(z_vec_long, x_vec, imag(cur_field_rep) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Full mode, w/ phase, (imag), mode %i of %i,\n k/b = %f, %i periods', i_mode, length(k), k(i_mode)/b, n_periods) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        subplot(1, 3, 2);
        imagesc(z_vec_long, x_vec, real(cur_field_rep) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Full mode, w/ phase, (real), mode %i of %i\n k/b = %f, %i periods', i_mode, length(k), k(i_mode)/b, n_periods) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        subplot(1, 3, 3);
        imagesc(z_vec_long, x_vec, abs(cur_field_rep) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Full mode, w/ phase, (abs), mode %i of %i,\n k/b = %f, %i periods', i_mode, length(k), k(i_mode)/b, n_periods) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
       
        % plot field, with phase, across several periods
        figure; set(gcf,'Position',[50 50 1700 900]);
        imagesc( imag(cur_field_rep) ); set(gca, 'Ydir','normal'); 
        title( sprintf('Full mode, w/ phase, (imag), mode %i of %i,\n k/b = %f, %i periods', i_mode, length(k), k(i_mode)/b, n_periods) );
        xlabel('z (nm)'); ylabel('x (nm)'); colorbar; colormap(jet);
        
        % calculate up/down power directivity
        i_top = 290; i_bot = 20;
        p_directivity = calculate_power_directivity( cur_field, dx, k(i_mode), i_top, i_bot );
        fprintf('Down/up power directivity over one unit cell is (Pdown/Pup) = %f\n', p_directivity );
        fprintf('Max down/up efficiency over one unit cell is  = %f percent\n', 100*(1/(1+1/p_directivity)) );
        
        % calculate up/down power directivity over many periods
        fprintf('Down/up power directivity over %i periods is (Pdown/Pup) = %f\n', n_periods, calculate_power_directivity( cur_field_rep, dx, k(i_mode), i_top, i_bot ) );
        
%         calculate_directivity( cur_field, dx, size(cur_field,1)-1, 1 )     % testing edge of domain
        
        % calculating angular directivity
        n_top = 1.4467; n_bot = 1;
        calc_angular_directivity( cur_field_rep, dx, i_top, i_bot, n_top, n_bot, lambda0 );
        
    end
%     pause(1);
%     close all;

    % calculate angle
    % k is kz
    theta_out = -acos( real(k)./k0 ) * 180/pi

end












