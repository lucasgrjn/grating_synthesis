function [ Ex, Ey ] = f_get_E_from_H_field( Hz, N, dx, dy, omega )
% returns the Ex and Ey components of TM 2d bloch mode
%
% currently ignores omega
%
% Inputs:
%   Hz
%       type: matrix
%       desc: H field, solved by modesolver. Includes phase.
%             Dimensions are Hz vs y vs x
%   N
%       type: matrix
%       desc: Index distribution. On same grid as Hz
%   dx
%   dy
%   omega

% coordinates
Nx          = size(N,2);
Ny          = size(N,1);
x           = 0:Nx-1;
y           = 0:Ny-1;
fx          = (-Nx/2:(Nx/2-1) )./Nx;                        % spatial freq vector along x
fy          = (-Ny/2:(Ny/2-1) )./Ny;                        % spatial freq vector along y
[ Fx, Fy ]  = meshgrid( fx, fy );

% taking spatial derivatives of Hz using fft
Hz_fft      = fftshift( fft2( ifftshift(Hz) ) );
dHz_dx_fft  = Hz_fft .* ( 1i * 2*pi*Fx );                       % x derivative, ft
dHz_dy_fft  = Hz_fft .* ( 1i * 2*pi*Fy );                       % y derivative, ft
dHz_dx      = fftshift( ifft2( ifftshift(dHz_dx_fft) ) );       % x deriv
dHz_dy      = fftshift( ifft2( ifftshift(dHz_dy_fft) ) );       % y deriv

% taking spatial derivatives using discrete derivative
dHz_dx = Hz(:,2:end) - Hz(:,1:end-1);
dHz_dy = Hz(2:end,:) - Hz(1:end-1,:);
% in this case, N must be placed on the half grid as well
x_halfgrid          = 0:1/2:Nx-1;
y_halfgrid          = 0:1/2:Ny-1;
[ X, Y ]            = meshgrid(x,y);
[ X_half, Y_half ]  = meshgrid( x_halfgrid, y_halfgrid );
N_halfgrid          = interp2( X, Y, N, X_half, Y_half );

% DEBUG plot N half grid
figure;
imagesc( N_halfgrid );
colorbar;
set(gca, 'ydir', 'normal');
title('DEBUG N half grid, linear interpolation');

% calc Ex and Ey (yee grid)
% Ex = dHz_dy./( 1i * N_halfgrid( 2:2:end-1, 1:2:end ).^2 );
% Ey = -dHz_dx./( 1i * N_halfgrid( 1:2:end, 2:2:end-1 ).^2 );
% calc Ex and Ey (yee grid)
Ex = dHz_dy;
Ey = -dHz_dx;

end

