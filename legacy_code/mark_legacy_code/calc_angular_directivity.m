function [ ] = calc_angular_directivity( field, dxz, i_top_row, i_bot_row, n_top, n_bot, lambda0 )
% calculates output angles
% 
% % TODO:
% right now does not spit out an output angle
% does plot the angular distribution for a slice through the top and a slice through the bottom
% what should I return from this function? 
% 
% 
% Inputs
%     field
%         : field distribution, including phase factor (exp(jkz) if bloch mode)
%           dimensions are (x,z)
%     dxz
%         : dicretization in x and z
%           units of nm
%     i_top_row
%         : index of top field row to slice and FFT
%     i_bot_row
%         : index of bottom field row to slice and FFT
%     n_top
%         : refractive index of top row
%     n_bot
%         : refractive index of bottom row
%     lambda0
%         : wavelength
%           units of nm

f0 = 1/lambda0; % units of spatial frequency

% get number of discrete points
[Nx, Nz] = size(field);

% create space and frequency vectors
x_vec   = dxz.*(0:Nx-1);             % space vector
z_vec   = dxz.*(0:Nz-1);             % space vector
fx_vec  = (-Nx/2:(Nx/2-1) ).*(1/(dxz*Nx));      % freq vector
fz_vec  = (-Nz/2:(Nz/2-1) ).*(1/(dxz*Nz));      % freq vector
angle_vec_bot               = (180/pi)*asin( fz_vec./(n_bot*f0) );      % draw a circ
real_a_indices_bot          = ( real(angle_vec_bot) == angle_vec_bot ); % use this to index the range of allowed free space modes
angle_vec_top               = (180/pi)*asin( fz_vec./(n_top*f0) );      % draw a circ
real_a_indices_top          = ( real(angle_vec_top) == angle_vec_top ); % use this to index the range of allowed free space modes

% FFT the bottom row
field_bot       = field( i_bot_row, : );
field_bot_fft   = fftshift( fft( ifftshift(field_bot) ) );   % essentially fz (spatial freq, = kz/2pi)
% FFT the top row
field_top       = field( i_top_row, : );
field_top_fft   = fftshift( fft( ifftshift(field_top) ) );   % essentially fz (spatial freq, = kz/2pi)

% DEBUG
% figure;
% plot( fz_vec, abs(field_bot_fft) );
% title('FFT of field (abs, k_z), bottom slice, vs fz');
% xlabel('spatial frequency f_z = k_z/2pi'); ylabel('FFT of field along z');
% makeFigureNice();

% DEBUG
figure;
plot( angle_vec_bot(real_a_indices_bot), abs(field_bot_fft(real_a_indices_bot)) );
title('FFT of field (abs, k_z), bottom slice, vs angle');
xlabel('\theta degrees'); ylabel('FFT of field along z');
makeFigureNice();

figure;
plot( fz_vec(real_a_indices_bot), abs(field_bot_fft(real_a_indices_bot)) );
title('FFT of field (abs, k_z), bottom slice, vs f_z, within the allowed free space modes');
xlabel('f_z'); ylabel('FFT of field along z');
makeFigureNice();

% DEBUG
figure;
plot( angle_vec_top(real_a_indices_top), abs(field_top_fft(real_a_indices_top)) );
title('FFT of field (abs, k_z), top slice, vs angle');
xlabel('\theta degrees'); ylabel('FFT of field along z');
makeFigureNice();

figure;
plot( fz_vec(real_a_indices_top), abs(field_top_fft(real_a_indices_top)) );
title('FFT of field (abs, k_z), top slice, vs f_z, within the allowed free space modes');
xlabel('f_z'); ylabel('FFT of field along z');
makeFigureNice();

% plotting both top and bottom
figure;
plot( angle_vec_bot(real_a_indices_bot), abs(field_bot_fft(real_a_indices_bot)) ); hold on;
plot( angle_vec_top(real_a_indices_top), abs(field_top_fft(real_a_indices_top)) );
title('FFT of field (abs, k_z) vs angle');
xlabel('\theta degrees'); ylabel('FFT of field along z');
legend('bottom', 'top');
makeFigureNice();

figure;
plot( fz_vec, 20*log10(abs(field_bot_fft)) ); hold on;
plot( fz_vec, 20*log10(abs(field_top_fft)) );
title('FFT of field (abs, k_z) vs f_z, dB scale');
xlabel('f_z'); ylabel('FFT of field along z');
legend('bottom', 'top');
makeFigureNice();

% plotting in polar
figure;
polarplot( angle_vec_bot(real_a_indices_bot)*pi/180, abs(field_bot_fft(real_a_indices_bot)) );
title('FFT of field (abs, k_z), bottom slice, vs angle');
set(gca, 'ThetaZeroLocation', 'bottom'); % rotate axes
% xlabel('\theta degrees'); ylabel('FFT of field along z');
% makeFigureNice();

figure;
polarplot( angle_vec_top(real_a_indices_top)*pi/180, abs(field_top_fft(real_a_indices_top)) );
title('FFT of field (abs, k_z), top slice, vs angle');
% xlabel('\theta degrees'); ylabel('FFT of field along z');
set(gca, 'ThetaZeroLocation', 'top'); % rotate axes
set(gca, 'ThetaDir', 'clockwise'); % rotate axes
% makeFigureNice();

% % how about along the tranverse direction?
% doesn't tell me much (not surprised...)
% i_col = 65;
% field_transv_slice = field( :, i_col); % just picking a random column
% field_transv_slice_fft = fftshift( fft( ifftshift( field_transv_slice ) ));
% 
% figure;
% plot( fx_vec, abs(field_transv_slice_fft) );
% title(sprintf('FFT of field (abs, k_x), transverse slice, vs fx, column = %i', i_col));
% xlabel('spatial frequency f_x = k_x/2pi'); ylabel('FFT of field along x');
% makeFigureNice();
% 
% figure;
% plot( angle_vec_fx(real_angle_indices_fx), abs(field_transv_slice_fft(real_angle_indices_fx)) );
% title(sprintf('FFT of field (abs, k_x), transverse slice, vs fx, column = %i', i_col));
% xlabel('\theta degrees'); ylabel('FFT of field along x');
% makeFigureNice();
% 
% % how about taking the FFT along the x dimension, then averaging the FFTs?
% field_transv_fft = fftshift( fft( ifftshift( field, 1 ) ), 1 ); % remember to fftshift the columns only
% field_transv_fft_avg = mean(field_transv_fft, 2);
% 
% figure;
% plot( fx_vec, abs(field_transv_fft_avg) );
% title( 'FFT of field (abs, k_x), all transverse slices averaged, vs fx' );
% xlabel('spatial frequency f_x = k_x/2pi'); ylabel('FFT of field along x');
% makeFigureNice();
% 
% figure;
% plot( angle_vec_fx(real_angle_indices_fx), abs(field_transv_fft_avg(real_angle_indices_fx)) );
% title( 'FFT of field (abs, k_x), all transverse slices averaged, vs fx' );
% xlabel('\theta degrees'); ylabel('FFT of field along x');
% makeFigureNice();

end

