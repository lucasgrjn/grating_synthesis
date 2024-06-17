function [synth_obj] = f_run_designspace_singlelevel( lambda, optimal_angle, disc )
% authors: bohan
% 
% Generates design space for example single level grating
%
% Inputs
%   lambda
%       wavelength (nm)
%   optimal_angle
%       desired output angle
%   disc
%       discretization (nm)

% initial settings
units               = 'nm';
index_clad          = 1.45;
y_domain_size       = 4000;
coupling_direction  = 'up';
data_notes          = ['lambda ' num2str(lambda) ' optimal angle ' num2str(optimal_angle)];

% display inputs
fprintf('Inputs are:\n');
fprintf('Wavelength: %f %s\n', lambda, units);
fprintf('Angle: %f degrees\n', optimal_angle);
fprintf('Discretization: %f nm\n', disc);

% handle to grating cell gen function
% one finicky detail is that the handle function must take the following arguments
%   dxy, background_index, y_domain_size, period, fill
% even if not all are used, or if you use other args instead
h_makeGratingCell = @( dxy, background_index, y_domain_size, ...
                         period, fill ) ...
                     f_makeGratingCell( dxy, ...
                                     period, fill );
                                         


% make synthesis object
synth_obj = c_synthGrating(   'discretization',    disc, ...
                              'units',             units,   ...
                              'lambda',            lambda, ...
                              'background_index',  index_clad,    ...
                              'coupling_index', index_clad, ...
                              'y_domain_size',     y_domain_size, ...
                              'optimal_angle',     optimal_angle, ...
                              'data_notes',        data_notes, ...
                              'coupling_direction', coupling_direction, ...
                              'h_makeGratingCell', h_makeGratingCell ...
                              );

% display Q for logging purposes
synth_obj

% design space fills to sweep
fills = 0.02:0.01:0.98;

% run design space generation
tic;
synth_obj = synth_obj.generate_design_space( fills );
toc;

fprintf('Design space generation sweep is done\n');

end


        