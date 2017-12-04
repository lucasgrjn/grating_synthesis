% Silicon electron and hole index change
% due to carrier concentration
%
% Syntax [dne,dnh,dae,dah] = soref_Si_carrier_to_index(deltaN,wavelength)
%
% Input:  deltaN - carrier concentration (cm^-3)
%         wavelength - one of selectable wavelengths ('1280' or '1550')
%
% Output: dne    - electron contribution to index change at given wavelength
%         dnh    - hole contribution to index change at given wavelength
%         dae    - electron contribution to absorption change at given wavelength
%         dah    - hole contribution to absorption change at given wavelength
%
% Reference: R.A. Soref and B.R. Bennett, IEEE J. Quantum Electron. QE-23, p. 123 (Jan 1987).

% Code updates:
% -------------
% 2011 Dec 22 - Added other wavelength data
% 2006 Nov 13 - First version.

function [dne,dnh,dae,dah] = soref_Si_carrier_to_index(deltaN,wavelength)
% Fit data from Soref plots (Fig. 11)
% From file: '20070114 - Milos - Data fits for Si free carrier data in Soref1987.ppt'
if strcmp(wavelength,'1280')
    nP1e = 1.066; nP2e = -22.51;
    nP1h = 0.8035; nP2h = -17.37;
    aP1e = 1.057; aP2e = -18.21;
    aP1h = 1.062; aP2h = -18.56;
elseif strcmp(wavelength,'1550')
    nP1e = 1.043; nP2e = -21.84;
    nP1h = 0.8067; nP2h = -17.23;
    aP1e = 1.098; aP2e = -18.71;
    aP1h = 1.105; aP2h = -19.17;
end

% carrierToIndex = @(DN,p1,p2) 10.^(p1*log10(DN)+p2);
carrierToIndex = @(DN,p1,p2) 10^p2 * DN.^p1;

dne = -carrierToIndex(deltaN,nP1e,nP2e);
dnh = -carrierToIndex(deltaN,nP1h,nP2h);
dae = 100*carrierToIndex(deltaN,aP1e,aP2e);
dah = 100*carrierToIndex(deltaN,aP1h,aP2h);
