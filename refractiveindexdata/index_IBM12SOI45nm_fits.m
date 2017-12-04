% Refractive index fits for IBM 12SOI 45nm process (used for EOS4 chip)
%
% n = index_IBM12SOI45nm_fits(lam_um, material)
%
% Outputs:
% n        - [vector] Refractive index vector
%
% Inputs:
% lam_um   - [vector] wavelength in microns
% material - [string] 'polySi' - polysilicon gate
%                     'STI' - shallow trench isolation under waveguide
%                     'PSG' - oxide over the waveguide and Si3N4 liner

function n = index_IBM12SOI45nm_fits(Lum, material)

switch(material)
    case 'polySi'
        n = 4.337 - 1.247 * Lum + 0.6795 * Lum.^2 - 0.1305 * Lum.^3;
    case 'STI'
        n = 1.472 - 0.02556 * Lum + 0.005342 * Lum.^2;
    case 'PSG'
        n = 1.488 - 0.02998 * Lum + 0.00676 * Lum.^2;
    otherwise
        error([mfilename ': No valid material specified.']);
end
