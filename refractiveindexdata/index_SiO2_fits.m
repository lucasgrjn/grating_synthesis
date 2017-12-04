% Silica (SiO2) index vs. wavelength (x-xum) and temperature (x-xK)
% M. Popovic, Apr 28, 2006
%
% Syntax:  [n nmin nmax] = index_SiO2_fits(lambda_um, temp_Kelvin)
%
% Input:    lambda_um      - [vector] wavelengths (accurate range 1.5-1.6um)
%           temp_Kelvin    - [vector] temperature (accurate range 85-920K)
%
% Output:   n              - nominal refractive index
%           nmin, nmax     - 95% confidence interval of my fit of
%           literature fit curves.
%
% Sources of data:
% ----------------
% G. Ghosh et al., "Temperature-dependent Sellmeier coefficients and
% chromatic dispersions for some optical fiber glasses," J. Lightwave
% Technol. vol. 12, no. 8, Aug. 1994, p.1338.

function [nn, nnmin, nnmax] = index_SiO2_fits(lam, temp)
if(nargin < 2) temp = 298.15; end   % Default room temperature
lam = lam(:); temp = temp(:);


TC = temp - 273.15;     % Celsius used in model
% Wavelength and temperature fit from above-cited paper
A = 6.90754e-6 * TC + 1.31552;
B = 2.35835e-5 * TC + 7.88404e-1;
C = 5.84758e-7 * TC + 1.10199e-2;
D = 5.48368e-7 * TC + 0.91316;
E = 100;
nn    = sqrt( A + B./(1 - C./lam.^2) + D./(1 - E./lam.^2) );  % Nominal n(lambda)
nnmin = NaN; %3.476 + (A-dA)*x + (B-dB)*x.^2;     % ..min and..
nnmax = NaN; %3.476 + (A+dA)*x + (B+dB)*x.^2;     % ..max of 95% conf interval
