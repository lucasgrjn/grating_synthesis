function [ thetas, mfds, coupling_vs_angle_mfd, overlap_vs_angle_mfd ] = f_coupling_vs_mfd_angle( Ez, Hx, lambda0, lambdas, n_background, center_angle, x, mfd, T  )
% Calculates coupling vs. angle and mfd and wavelength
% Same as the script but in functional format
%
% Inputs:
%   Ez
%       type: double, matrix
%       desc: E field (z component) at bottom monitor
%             dimensions are x vs. wavelength
%   Hx
%       type: double, matrix
%       desc: H field (x component) at bottom monitor
%             dimensions are x vs. wavelength
%   lambda0
%       type: double, scalar
%       desc: Desired wavelength to compute coupling at
%   lambdas
%       type: double, array
%       desc: simulated wavelengths
%   n_background
%       type: double, scalar
%       desc: index of refraction of cladding at monitor plane
%   center_angle
%       type: double, scalar
%       desc: coupling angle which the grating was designed for
%   x
%       type: double, array
%       desc: x coordinates, in meters
%   mfd
%       type: double, scalar
%       desc: mode field diameter, in meters
%   T
%       type: double, array
%       desc: transmission from bottom field monitor vs. wavelength
%
% Outputs:
%   thetas
%   mfds
%   coupling_vs_angle_mfd
%   overlap_vs_angle_mfd

% squeeze data
Ez = squeeze( Ez );                                               % dimensions x vs. freq
Hx = squeeze( Hx );                                               % dimensions x vs. freq

% % grab center wavelength (or closest to it)
[~, indx_center_wl]   = min( abs( lambdas - lambda0 ) ); 


% grab field @ center wavelength
Ez_center_wl    = Ez( :, indx_center_wl );
Hx_center_wl    = Hx( :, indx_center_wl );       % currently unused in the overlap formula
T_center_wl     = T( indx_center_wl );

% mfds
mfds = linspace( mfd*0.5, mfd*1.5, 20 );

% calculate coupling vs. angle
d0              = 0;
unit_scaling    = 1.0;
nclad           = n_background;
w0              = mfd/2; % 10.4e-6 / 2;                                              % MFD/2
x_fiber         = x - x( round(end/2) );
thetas          = linspace( center_angle - 25, center_angle + 25, 101 );
% thetas          = -15;

% saving variables
coupling_vs_angle_mfd    = zeros( length(mfds), length(thetas) );          % dimensions coupling vs. mfd vs. theta
overlap_vs_angle_mfd     = zeros( length(mfds), length(thetas) );          % dimensions overlap vs. mfd vs. theta

% freqs = freqs(end/2);
for i_theta = 1:length(thetas)
 
    % calculate overlap vs MFD
    field_overlap = zeros( size(Ez,1), length(mfds) );                              % dimensions x vs. mfd
    for i_mfd = 1:length(mfds)
        
        % make fiber mode
        w0                   = mfds(i_mfd)/2;
        [Ez_fiber, Hx_fiber] = f_fiberModeGaussian( w0, lambdas(indx_center_wl), x_fiber, unit_scaling, -thetas(i_theta), d0, nclad );
        
        % calc overlap
        field_overlap(:,i_mfd) = f_field_position_overlap_fullvec_1D( ...
                                            struct( 'z', Ez_center_wl ), struct( 'x', Hx_center_wl ), ...
                                            struct( 'z', Ez_fiber ), struct( 'x', Hx_fiber ), ...
                                            'y');
        
    end
    
    coupling_vs_angle_mfd(:,i_theta) = T_center_wl .* max( abs( field_overlap ) ).';
    overlap_vs_angle_mfd(:,i_theta)  = max( abs( field_overlap ) ).';
    
end

end