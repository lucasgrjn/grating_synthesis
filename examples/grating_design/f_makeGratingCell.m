function GC = f_makeGratingCell( ...
    dxy, period, dutycycle )
% makes and returns a c_twoLevelGratingCell object
% with "generic" SOI process
%
% units are in nm, specifically for this method
% 
% inputs:
%   dxy
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y, in units of 'units'
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
%   period
%       type: double, scalar
%       desc: period of the grating cell, in units defined by synth_obj.units
%   dutycycle
%       type: double, scalar
%       desc: ratio of waveguide to period
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object

% material indices
n_sio2 = 1.45;
n_si = 3.47;

% waveguide thicknesses
wg_thick = 220;
paretch_depth = 110;

% set domain 
y_domain_size = 4e3;
domain_size = [ y_domain_size, period ];

% make single level grating cell
GC = c_gratingCell( 'discretization', dxy, ...
                    'domain_size', domain_size, ...
                    'background_index', n_sio2, ...
                    'numcells', 10 );
     
% get air thickness
domain_y_half = round( (domain_size(1)/2) /GC.dy) * GC.dy;

% layer min heights
wg_miny = domain_y_half - wg_thick/2;

% draw the geometry
% in this case, a 220 nm thick waveguide with a 110 nm partial etch

% first add silicon waveguide layer
GC = GC.addLayer(wg_miny, wg_thick, n_si);

% then add etched region
etch_miny = wg_miny + wg_thick - paretch_depth;
etch_width = (1 - dutycycle) * period;
GC = GC.addRect(0, etch_miny, etch_width, paretch_depth, n_sio2);

% set top and bottom bounds for waveguide region
% this is important for choosing the right mode
GC.wg_min_y = wg_miny; % bottom position of wg
GC.wg_max_y = wg_miny + wg_thick; % top position of wg

             
end




