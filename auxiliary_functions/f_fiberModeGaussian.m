function [Ez, Hx, E_tm, H_tm, mfd] = f_fiberModeGaussian( w0, lambda0, xvec, unit_scale, theta, d0, nclad )
% somewhat adapted from Cale's code
%
% Generate Gaussian-beam mode profile at a plane through y = d0 at angle of
% theta (in degrees)
%
% NOTE that the angle theta is with respect to the positive y axis. When
% overlapping downwards, a "positive theta" in our typical definition (from
% the downwards normal) is actually a negative theta with respect to the
% positive y axis.
%
% now solves for TM too, but is kind of hacky
%
% (using H.A. Haus, Waves & Fields, Chapter 5)
% CMG November 21st, 2014
%
%
% Inputs:   
%   w0  
%       type: double, scalar
%       desc: 1/e beam RADIUS at waist, in units 'units'
%           
%   lambda0
%       type: double, scalar
%       desc: free space wavelength, in units 'units'
%
%   xvec
%       type: double, array
%       desc: coordinates along direction of propagation of
%             grating, in units 'units'
%
%   unit_scale
%       type: double, scalar
%       desc: scaling factor, such that 'units' * unit_scale = meters
%             I don't think this is even necessary tho!
%               The only dependence is calculating omega0, which should be
%               spatial unit free anyways... IDK i'll figure this out later
%           
%   theta 
%       type: double, scalar
%       desc: angle from normal in degrees
%
%   d0  
%       type: double, scalar
%       desc: distance from beam waist to slice
%
%   nclad 
%       type: double, scalar
%       desc: cladding index
%   pol
%       type: str
%       desc: either 'TE' or 'TM'
%             optional for now, and defaults to TE
%
%
% Outputs: 
%   u
%       type: double, array
%       desc: returned slice of gaussian beam, normalized to total
%             power
%   E_tm and H_tm
%       this is a new addition - saves E and H fields for "TM" fiber
%

% % default polarization to TE
% if nargin < 8
%     pol = 'TE';
% end


% Constants, in units of meters
c       = 3e8;                                  % m/s
mu0     = 4*pi * 1e-7;                          % H/m
omega0  = 2*pi*c/( lambda0 * unit_scale);       % rad/s.
lambda  = lambda0 * unit_scale / nclad;         % wavelength in cladding, units m
k       = 2*pi/lambda;                          % 1/m
w0      = w0 * unit_scale;                      % [meters] radius
d0      = d0 * unit_scale;                      % [meters] offset

% Convert to radians
theta = (pi/180)*theta;

% Scale coordinates
xvec = xvec * unit_scale;                                              % units m

% coordinates in fiber frame
xprime = xvec.*cos(theta) - d0*sin(theta);
zprime = xvec.*sin(theta) + d0*cos(theta);

% b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)
b = k*w0^2/2;                                                                                   

% Equation (5.2) in Haus [1/meters]
u00 =   1j .* sqrt(k*b/pi) .* ( 1./(zprime + 1j.*b) ).*...
    exp( -1j.*k.*( xprime.^2 )./( 2*(zprime + 1j.*b) ) );     

% normalize the mode to intensity, makes things nicer
dx      = xvec(2) - xvec(1);
u00     = u00/sqrt( dx * sum( abs( u00 ).^2 ) );


% calculate fields
E_yprime = -1i * omega0 * u00 .* exp( -1i * k * zprime );
H_xprime = ( -1i * k/mu0 ) * u00 .* exp( -1i * k * zprime );
H_zprime = ( -xprime./( zprime + 1i*b ) ) .* H_xprime; 

% project fields back into original coordinates
Ez = E_yprime;                                                              % this one is used (TE pol grating)
Hx = H_xprime * cos(theta) + H_zprime * sin(theta) ;

% if TM, project fields differently
% pretty damn confusing, but check 2019 10 - clo / 2019 10 07 - tm grating
% for derivation
E_tm.x = E_yprime * cos(theta);
E_tm.y = E_yprime * sin(theta); 
H_tm.z = H_xprime;
H_tm.x = H_zprime * sin(theta);
H_tm.y = H_zprime * cos(theta);

% TEMP using formula to calc MFD
r       = xprime(xprime >= 0);
E_temp  = E_yprime(xprime >= 0);
dr      = r(2:end)-r(1:end-1);      % calc dr for integral
r       = r(1:end-1);               % put on dr grid
E_temp  = E_temp(1:end-1);          % put on dr grid

% calc numerator
numer = sum( (E_temp.^2) .* (r.^3) .* dr );

% calc denom
denom = sum( (E_temp.^2) .* r .* dr );

% calc mfd
mfd = 2*sqrt(2*(numer/denom));



end     % end fiberModeGaussian()












































