% unidirectional grating
% trying to recreate jelena's results
% now trying to recreate marks' results

clear; clc; close all;
tic;
% all units nm

addpath('refractiveindexdata');

% ------------------------------------------------------------------------|
% Parameters
% ------------------------------------------------------------------------|

lambda0     = 1310;     % in nm
k0          = 2*pi/lambda0;
temp_kelvin = 273;  % room temperature

% indices, 
[n_cSi, ~, ~]   = index_Si_fits(lambda0/1000, temp_kelvin);
n_pSi           = index_IBM12SOI45nm_fits(lambda0/1000, 'polySi');
n_SiN           = 2.0;
[n_SiO2, ~, ~]  = index_SiO2_fits(lambda0/1000, temp_kelvin);

% dimensions
period  = 900; % 650;       % grating period
w_t     = 500;       % width top of grating
w_b     = 0;       % width bottom of grating
t       = 80;       % thickness of grating teeth and SiN layer
w_o     = 150;
t_air   = 0;    % 700;
t_SiO2_bot = 1000;   % thickness of SiO2 layer bottom
t_SiO2_top = 1000;   % thickness of SiO2 layer top

% discretization
dxz = 10;

% ------------------------------------------------------------------------|
% Make grating index
% ------------------------------------------------------------------------|

[ n, x_vec, z_vec ] = make_grating_cell( period, dxz, [ n_cSi, n_pSi, n_SiN, n_SiO2 ],...
                                    w_t, w_b, w_o, t, t_air, t_SiO2_bot, t_SiO2_top );

% display
figure; imagesc(z_vec, x_vec, n); set(gca,'Ydir','normal'); colorbar; axis equal;
colormap(jet); title('Index of unit cell');
xlabel('z (nm)'); ylabel('x (nm)');

% ------------------------------------------------------------------------|
% Run modesolver
% ------------------------------------------------------------------------|

% inputs
modes   = 10;
BC      = 0;         % 0 is pec, 1 is pmc, doesn't matter if using PML's tho

% pml options
h_pml       = 500; % nm
pml_str     = 600;
pml_order   = 2;
PML_options = [1, h_pml, pml_str, pml_order];

n_tries = 1;

% try a bunch of guess k's
err_flag = false;
% theta_try = (20) * pi/180;
guessk = pi/period;
% guessk = 0.07362*2*pi/Lambda + 0.0003619i;
% guessk = 0
% guessk = 0.0008 + 2*pi/Lambda;
% guessk = k0*1.1
% guessk = 7.4896e-04
% guessk = k0*cos(theta_try);
% guessk = 0.0026 - 0.0003i;
% guessk = 0.0012 + 0.0004i;
for ii = 1:n_tries
    
    fprintf('try %i of %i\n',ii,n_tries);
    
    try
        % run modesolver
        [Phi_1D, k] = complexk_mode_solver_2D_PML(n,dxz,k0,modes,guessk,BC,PML_options);
    catch
%         guessk = k0*(1 + 0.1*ii);
%         continue;
    end
    
    break;
    
end

% ------------------------------------------------------------------------|
% show results
% ------------------------------------------------------------------------|

% plot fields
figure_pos = [50 50 1300 600];

for ii = 1:length(k)
    
    % reshape field
    field = reshape(Phi_1D(:,ii), fliplr(size(n))); % dimensions (z, x) where z = direction of propagation
    field = field.';
    field = field./max(abs(field(:)));
    field = fliplr( field );    % i think the field is z flipped
    
    
    % add in the phase term
    k_phase         = exp(1i*k(ii).*z_vec);
    field_w_phase   = field.*repmat(k_phase,length(x_vec),1);
%     figure;
%     plot(

    % plot field, no phase ----------------------
    % imaginary
    figure; 
    set(gcf,'Position', figure_pos);
    subplot(1,3,1);  
    imagesc(z_vec, x_vec, imag(field)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, envelope only, imag', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % real
    subplot(1,3,2); 
    imagesc(z_vec, x_vec, real(field)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, envelope only, real', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % abs
    subplot(1,3,3); 
    imagesc(z_vec, x_vec, abs(field)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, envelope only, abs', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % plot field, w/ phase ----------------------
    % imag
    figure;  
    set(gcf,'Position', figure_pos);
    subplot(1,3,1); 
    imagesc(z_vec, x_vec, imag(field_w_phase)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, imag', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % real
    subplot(1,3,2); 
    imagesc(z_vec, x_vec, real(field_w_phase)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, real', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % abs
    subplot(1,3,3); 
    imagesc(z_vec, x_vec, abs(field_w_phase)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, abs', ii, length(k)) );
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % plot field with phase multiple periods ----------------------
    % extend field
    n_periods           = 10;
    field               = repmat(field,1,n_periods);
    z_vec_long          = linspace(0, n_periods*period, size(field,2) );
    k_phase             = exp(1i*k(ii).*z_vec_long);
    field_w_phase_long  = field.*repmat(k_phase,length(x_vec),1);
    
    % imag
    figure;  set(gcf,'Position', figure_pos);
    subplot(1,3,1); 
    imagesc(z_vec_long, x_vec, imag(field_w_phase_long)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, %i periods, imag', ii, length(k), n_periods ));
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % real
    subplot(1,3,2); 
    imagesc(z_vec_long, x_vec, real(field_w_phase_long)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, %i periods, real', ii, length(k), n_periods ));
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );
    
    % abs
    subplot(1,3,3); 
    imagesc(z_vec_long, x_vec, abs(field_w_phase_long)); 
    set(gca, 'Ydir','normal'); colorbar;
    title( sprintf('Mode: %i of %i, with plane wave phase, %i periods, abs', ii, length(k), n_periods ));
    caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]);
    xlabel('z (nm)'); ylabel('x (nm)');
    colormap( redbluehilight() );


end

% display wavevector value
k*period/(2*pi)

% k/k0
% theta_out = (180/pi)*atan( real(k)./real(sqrt( (k0*n_cSi)^2 - k.^2)) )

toc









