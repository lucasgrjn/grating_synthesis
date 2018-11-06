function GC = f_makeGratingCell_basic( dxy, units, lambda, background_index, y_domain_size, ...
                                         period, fill_top, fill_bot, offset_ratio )
% makes and returns a c_twoLevelGratingCell object
% Basic layer stack, 2 100nm thick Si layers that can be fully etched
% Fully cladded by SiO2
% 
% inputs:
%   dxy
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y, in units of 'units'
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
%   units
%       type: string
%       desc: name and scaling of spatial units, supports 'm'
%             (meters), 'mm' (millimeters), 'um' (microns), 'nm'
%             (nanometers)
%   lambda
%       type: double, scalar
%       desc: wavelength to solve at, in units 'units'
%   background_index
%       type: double, scalar
%       desc: value of background index, actually unused here
%   y_domain_size
%       type: double, scalar
%       desc: vertical/transverse domain size
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

% defaults
background_index = 1.45;

% set domain 
domain_size = [ y_domain_size, period ];

% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );

% make 2 level grating cell
GC = c_twoLevelGratingCell( 'discretization', dxy, ...
                            'units', units, ...
                            'lambda', lambda, ...
                            'domain_size', domain_size, ...
                            'background_index', background_index, ...
                            'numcells', 10 );

% index of refraction
n_SiO2  = 1.45;
n_Si    = 3.47;

% define layer thicknesses (nm)
domain_y_half   = round( (domain_size(1)/2) /GC.dy) * GC.dy;
t_Si_bottom     = 100;
t_Si_top        = 100;
                        
% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = [ t_Si_top, t_Si_bottom ];
wg_min_y        = [ domain_y_half-t_Si_bottom, domain_y_half ];
wgs_duty_cycles = [ fill_top, fill_bot ];
wgs_offsets     = [ 0, offset_ratio*period ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, [ n_Si, n_Si ], ...
                                        wgs_duty_cycles, wgs_offsets );
             
end


% -------------------------------------------------------------------------
% Auxiliary functions
% -------------------------------------------------------------------------



