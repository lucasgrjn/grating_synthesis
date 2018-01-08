function [ FOM ] = f_merit_period_offset( inputs, synthGrating_obj, angle, weights, fill_factors )
% authors: bohan zhang
% 
% Merit function that will be used to optimize grating cell's period and
% offset
%
% currently assuming we're propagating downwards
% 
% inputs:
%   inputs
%       type: 1x2 array
%       desc: [period, offset], where period is in units of nm/1000
%   synthGrating_obj
%   angle - angle in deg
%   weights - 1x2 array to weigh the two objectives

% parse inputs
period      = inputs(1)*1000;            % must convert to nm
offset      = inputs(2);
fill_top    = fill_factors(1);
fill_bot    = fill_factors(2);

% make grating coupler object
GC = synthGrating_obj.h_makeGratingCell( synthGrating_obj.convertObjToStruct(), period, fill_top, fill_bot, offset );

% simulation settings
num_modes   = 5;
BC          = 0;     % 0 for PEC, 1 for PMC
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
pml_options = [ 1, 200, 20, 2 ];

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options );

% minimize the FOM
FOM =   weights(1) * abs( angle - GC.max_angle_down )/angle + ...
        weights(2) * log10(GC.directivity);

end

