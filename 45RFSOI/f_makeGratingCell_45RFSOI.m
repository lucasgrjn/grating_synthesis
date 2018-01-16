function GC = f_makeGratingCell_45RFSOI( synth_obj, period, fill_top, fill_bot, offset_ratio )
% makes and returns a c_twoLevelGratingCell object
% with the 45RFSOI process parameters
%
% units are in nm, specifically for this method
% 
% inputs:
%   synth_obj
%       type: c_synthGrating object AS STRUCT
%       desc: c_synthGrating object AS STRUCT
%   period
%       type: double, scalar
%       desc: period of the grating cell, in units defined by synth_obj.units
%   fill_top
%       type: double, scalar
%       desc: ratio of top layer to period
%   fill_bot
%       type: double, scalar
%       desc: ratio of bottom layer to bottom layer
%   offset_ratio
%       type: double, scalar
%       desc: ratio of bottom layer offset to period
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object
%
% example:
%   period          = 800;
%   fill_top        = 0.8;
%   fill_bot        = 0.6;
%   offset_ratio    = 0.0;
%   GC              = f_makeGratingCell_45RFSOI( Q.convertObjToStruct(), period, fill_top, fill_bot, offset_ratio );
%       Then you can plot it and look at what the dielectric looks like:
%   GC.plotIndex()


% set domain 
domain_size     = synth_obj.domain_size;
domain_size(2)  = period;

% make grating cell
GC = c_twoLevelGratingCell( 'discretization', synth_obj.discretization, ...
                            'units', synth_obj.units.name, ...
                            'lambda', synth_obj.lambda, ...
                            'domain_size', domain_size, ...
                            'background_index', synth_obj.background_index, ...
                            'numcells', 10 );

% index of refraction
lambda_um = synth_obj.lambda * synth_obj.units.scale * 1e6;
n_SiO2  = index_SiO2_fits(lambda_um);
n_SiN   = index_SiN(lambda_um);
n_cSi   = index_Si_fits(lambda_um);
n_pSi   = index_IBM12SOI45nm_fits(lambda_um, 'polySi');
                        

% grab pml parameters
% [ yes/no, length in nm, strength, pml poly order ]
pml_opts    = synth_obj.modesolver_opts.pml_options; 
pml_length  = pml_opts(2);

% define layer thicknesses
t_air       = 500;
t_SiO2_bot  = 500;
t_SiN       = 70;
t_cSi       = 80;
t_pSi       = 80;

% draw layers
GC = GC.addLayer( t_air, domain_size(1)-t_SiO2_bot, n_SiO2 );     % add in SiO2
                        
% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = [ t_pSi, t_cSi ];
wg_min_y        = [ t_air+t_SiO2_bot+wg_thick(2), t_air+t_SiO2_bot ];
% wgs_duty_cycles = [ fill*ratio, fill ];
wgs_duty_cycles = [ fill_top, fill_bot ];
wgs_offsets     = [ 0, offset_ratio*period ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, [ n_pSi, n_cSi ], ...
                                        wgs_duty_cycles, wgs_offsets );


% draw SiN layer
top_layer_length    = fill_top*period;
% top of Si
GC = GC.addRect( 0, ...
                 t_air + t_SiO2_bot + wg_thick(2) + wg_thick(1), ...
                 top_layer_length + t_SiN, ...
                 t_SiN, ...
                 n_SiN );
% sidewall 1
GC = GC.addRect( top_layer_length, ...
                 t_air + t_SiO2_bot + wg_thick(2), ...
                 t_SiN, ...
                 wg_thick(1), ...
                 n_SiN );
% top of SiO2
GC = GC.addRect( top_layer_length, ...
                 t_air + t_SiO2_bot + wg_thick(2), ...
                 period - top_layer_length, ...
                 t_SiN, ...
                 n_SiN );
% sidewall2
if (top_layer_length <= period - t_SiN)
    % only draw 2nd sidewall in scenario that there is enough space for it
    GC = GC.addRect( period - t_SiN, ...
                     t_air + t_SiO2_bot + wg_thick(2), ...
                     t_SiN, ...
                     t_SiN + wg_thick(1), ...
                     n_SiN );
end
             
end


% -------------------------------------------------------------------------
% Auxiliary functions
% -------------------------------------------------------------------------



function n = index_IBM12SOI45nm_fits(Lum, material)
% Refractive index fits for IBM 12SOI 45nm process (used for EOS4 chip)
%
% author: unknown
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

end         % end index_IBM12SOI45nm_fits()



function [nSi, nSimin, nSimax] = index_Si_fits(lam, temp)
% Silicon index vs. wavelength (1.5-1.6um) and temperature (85-920K)
% (simple fits of more complicated fits in source papers)
% author: M. Popovic, Apr 1, 2006
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

end         % end index_Si_fits()



function [nn, nnmin, nnmax] = index_SiO2_fits(lam, temp)
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

end         % end index_SiO2_fits()




function nn = index_SiN(lam)
% SiN refractive index fits to expt'al data
% http://www.filmetrics.com/refractive-index-database/Si3N4/Silicon-Nitride-SiN
%
% Syntax:  nn = index_SiN(lambda_um)
%
% Input:    lambda_um   [1-vector] wavelengths in microns; if blank, a plot is given.
%          
% Output:   nn          [1-vector] refractive index
%
% Cale Gentry, June 24, 2014

p = [-0.007287569871325   0.036401331486906  -0.079722168861525   2.052602158358897];

nn = p(1)*lam.^3+p(2)*lam.^2+p(3)*lam+p(4);

end





