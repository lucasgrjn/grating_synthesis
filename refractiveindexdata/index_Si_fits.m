% Silicon index vs. wavelength (1.5-1.6um) and temperature (85-920K)
% (simple fits of more complicated fits in source papers)
% M. Popovic, Apr 1, 2006
%
% Syntax:  [nSi nSimin nSimax] = index_Si_fits(lambda_um, temp_Kelvin)
%
% Input:    lambda_um      - [vector] wavelengths (accurate range 1.5-1.6um)
%           temp_Kelvin    - [vector] temperature (accurate range 85-920K)
%
% Output:   nSi            - nominal refractive index
%           nSimin, nSimax - 95% confidence interval of my fit of
%           literature fit curves.
%
% Sources of data:
% ----------------
% D.F. Edwards, "Silicon (Si)" in Palik, Handbook of Optical Consts, 1985.
%
% J.A. McCaulley et al., Phys. Rev. B 49, p7408 (Mar 1994): ‘Temperature
% dependence of the near-IR refractive index of Silicon, GaAs and InP’.

function [nSi, nSimin, nSimax] = index_Si_fits(lam, temp)
if(nargin < 2) temp = 298.15; end   % Default room temperature
lam = lam(:); temp = temp(:);

% Wavelength
A = -0.07805; dA = 0.0005;
B = 0.082; dB = 0.002;
x = lam - 1.55;
nn    = 3.476 + A*x + B*x.^2;               % Nominal n(lambda)
nnmin = 3.476 + (A-dA)*x + (B-dB)*x.^2;     % ..min and..
nnmax = 3.476 + (A+dA)*x + (B+dB)*x.^2;     % ..max of 95% conf interval

% Temperature
A = 0.01587; dA = 0.00001;
B = 0.002392; dB = 0.000005;
To = 298.15;  x = temp/To - 1;
nTovnTo = 1 + A*x + B*x.^2;
nTovnTomin = 1 + (A-dA)*x + (B-dB)*x.^2;
nTovnTomax = 1 + (A+dA)*x + (B+dB)*x.^2;

% Combine wavelength and temp, and set output
nSi = nn .* nTovnTo;  nSimin = nnmin .* nTovnTomin;  nSimax = nnmax .* nTovnTomax;
