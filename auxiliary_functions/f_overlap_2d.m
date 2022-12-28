function [overlap_vs_wl_at_best_pos, indx_max_row, indx_max_col] = ...
        f_overlap_2d(Ex, Ey, Hx, Hy, w0, lambda, theta, x, y, nclad, coupling_direction, pol )
% helper function, calculates overlap in 2D at one angle, picks out the wavelength that
% has the best overlap, and returns the overlap vs wavelength when the fiber is
% positioned for coupling at that wavelength
%
% only for coupling to gaussian
% supports TE or TM coupling
%
% Inputs:
%   Ex
%       TE field, dimensions are x vs y vs frequency
%   Ey
%       TM field dimensions are x vs y vs frequency
%   Hx, Hy, similar to Ex and Ey
%   w0
%       MFD/2
%   lambda
%       wavelength vector
%   theta
%       angle to overlap at, scalar
%   x, y
%       positional vectors
%   nclad
%       cladding index, scalar
%   coupling_direction
%       either 'up' or 'down'
%   pol
%       either 'te' or 'tm'
%
% the input fields come from fdtd after they've been loaded, see this
% script for example: 
% D:\Google Drive\research\popovic group\projects\polarization insensitive grating\code\2019 10 18 - 45CLO TE TM degenerate grating\s_load_3dfdtd_te_tm_grating_updated_211106.m

%   
d0              = 0;
x_fiber         = x - x(round(end/2));
y_fiber         = y - y(round(end/2));

% calculate overlap using function
field_overlap       = zeros( size(Ex) );                              % dimensions x vs. freq
max_overlap         = zeros( size(lambda) );
for ii = 1:length(lambda)
    
    fprintf('freq loop %i of %i\n', ii, length(lambda) );

    % make fiber mode
    [E_fiber, H_fiber]  = f_fiberModeGaussian_2D( w0, lambda(ii), x_fiber, y_fiber, -theta, d0, nclad );
    [E_fiber, H_fiber]  = f_shape_gaussian_direction_pol( E_fiber, H_fiber, coupling_direction, pol );
        
    field_overlap(:,:,ii)     = f_field_position_overlap_fullvec_2D( ...
                                            struct( 'x', Ex(:,:,ii), 'y', Ey(:,:,ii) ), ...
                                            struct( 'x', Hx(:,:,ii), 'y', Hy(:,:,ii) ), ...
                                            E_fiber, H_fiber, ...
                                            'z');
    max_overlap(ii) = max( max( abs( field_overlap(:,:,ii) ) ) );
    
end

% pick out the coupling eff at the position for best coupling
[~, indx_wl_best]   = max(max_overlap);
[ indx_max_row, indx_max_col ] = ...
    find( max( max( abs(field_overlap(:,:,indx_wl_best)) ) ) == abs(field_overlap(:,:,indx_wl_best)) );
overlap_vs_wl_at_best_pos = field_overlap( indx_max_row, indx_max_col, : );

end

