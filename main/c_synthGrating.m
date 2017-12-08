classdef c_synthGrating
% Synthesizes a 2-level grating in an arbitrary process
%
% Authors: bohan zhang
%
%
% Based on Mark/Jelena's synthesis suite/pipeline
%
%
% Prerequisites/dependencies
%   - c_twoLevelGratingCell.m
%   - the utility folder
%
% Running the parameter sweep:
%   The user is allowed to define their own custom grating unit cell
%   drawing function.
%   HOWEVER, this function MUST have the following inputs and outputs, IN
%   ORDER:
%       function GC = your_makeGratingCell_function( synth_obj, period, fill, ratio, offset_ratio )
%           % makes and returns a c_twoLevelGratingCell object
% 
%           inputs:
%               synth_obj
%               type: c_synthGrating object AS STRUCT
%               desc: c_synthGrating object AS STRUCT
%           period
%               type: double, scalar
%               desc: period of the grating cell
%           fill
%               type: double, scalar
%               desc: ratio of bottom layer to period
%           ratio
%               type: double, scalar
%               desc: ratio of top layer to bottom layer
%           offset_ratio
%               type: double, scalar
%               desc: ratio of bottom layer offset to period
%
%   As an example:
%       see the default makeGratingCell() function that's included at the
%       end of this class definition. It makes a simple waveguide + air
%       cladding grating, no fancyness
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object
%
%
% Inputs to constructor:
%   Inputs are name-value pairs:
%   'discretization'
%       type: double, scalar
%       desc: discretization along x and y, in units of 'units'
%
%   'units'
%       type: string
%       desc: name and scaling of spatial units, supports 'm'
%             (meters), 'mm' (millimeters), 'um' (microns), 'nm'
%             (nanometers)
%
%   'lambda'
%       type: double, scalar
%       desc: wavelength to solve at, in units 'units'
%
%   'background_index'
%       type: double, scalar
%       desc: value of background index
%
%   'domain_size'
%       type: 1x2 array, double
%       desc: domain size, [ y height, x length ]
%
%   'period_vec'
%       type: double, array
%       desc: vector of periods to sweep, in units 'units'
%
%   'offset_vec'
%       type: double, array
%       desc: vector of offsets to sweep, as ratio of offset length/period
%
%   'ratio_vec'
%       type: double, array
%       desc: vector of ratio of bottom tooth to top tooth = bot tooth/top
%             tooth
%
%   'fill_vec'
%       type: double, array
%       desc: vector of ratio of top tooth to period = top tooth/period
%
%   'optimal_angle'
%       type: double, scalar
%       desc: desired output angle, in deg
%
%   'waveguide_index'
%       type: double, 1x2 array
%       desc: waveguide indexes [ index of top tooth, index of bottom
%             tooth ]
%
%   'waveguide_thicks'
%       type: double, 1x2 array
%       desc: waveguide thicknesses [ top thickness, bottom thickness ] 
%
%   'coupling_direction'
%       type: string
%       desc: direction of output light, 'up' or 'down'
%
%   'data_directory'
%       type: string
%       desc: path to data save directory
%
%   'data_filename'
%       type: string
%       desc: name of data file to save to/load from
%
%   'data_notes'
%       type: string
%       desc: optional verbose notes/descriptor for this simulation
%
%   'data_mode'
%       type: string
%       desc: flag to set data loading mode.
%             use 'new' to start a fresh simulation from scratch or 'load'
%             to load previously simulated data
%
%   'num_par_workers'
%       type: int, scalar
%       desc: number of parallel workers to use when running sweep

    properties
        % OLD properties
%         fillVector
%         ratioVector
%         periodVector
%         offsetVector
%         periodMatrix %note: these are really 4 dimensional arrays. need to change naming from Matrix to Array but will break some things
%         offsetMatrix
%         scatteringStrengthMatrix
%         directivityMatrix
%         angleMatrix
%         lambda0
%         P
%         inputs

        discretization;     % dx and dy
        units;              % units, verbose, 'm' or 'mm', or 'um', or 'nm'
                            % has fields 'name' and 'scale'
        lambda;             % center wavelength
        background_index;   % background index
        domain_size;        % domain size, [ y size (height), x size (length) ]
        period_vec;         % periods to sweep
        offset_vec;         % offsets to sweep
%         ratio_vec;          % ratios of top to bottom teeth lengths to sweep
%         fill_vec;           % fill ratios of top teeth to sweep
        fill_top_vec;       % top tooth fill ratio
        fill_bot_vec;       % bot tooth fill ratio
        optimal_angle;      % angle to optimize for, deviation from the normal, in deg.
        inputs;             % saves input settings for user reference

        start_time;         % time when object was created, 'YEAR-month-day hour-min-sec'
        
        % i'm thinking of getting rid of these two properties, since
        % they'll be defined by the end-user's custom grating cell function
        % anyways
        waveguide_index;    % [ <index of top tooth>, <index of bottom tooth> ]
        waveguide_thicks;   % [ <thickness of top tooth>, <thickness of bottom tooth> ]
        
        coupling_direction; % direction of coupling, either 'up', or 'down'
                            % defaults to 'down'
                            
        data_directory;     % path of data directory
        data_filename;      % name of data file
        data_notes;         % verbose notes of what current sweep is doing
        data_mode;          % flag telling whether grating is run from scratch or run from previous data
                            % either 'new' or 'load'
        
        sweep_results;      % struct holding results of parameter sweep
                            % tensors have dimensions ( fill, ratio, period, offset )
                            % AS OF 2017/12/08 CHANGED TO BE (fill top,
                            % fill bot, period, offset)
                            % fields are: fill_tensor, ratio_tensor, offset_tensor, period_tensor, 
                            % scatter_strengths, directivities, angles
        
        u;                  % gaussian profile, not sure if this will stay a property tho
        
        num_par_workers;    % number of parallel workers to use
        
        modesolver_opts;    % STRUCT that stores the modesolver options
                            % CURRENTLY hardcoded.
                            
        final_design;       % STRUCT that stores the final design parameters
                            % as well as performance values
                            
        h_makeGratingCell;  % handle to the grating cell making function
                            
    end
    
    methods
        
        function obj = c_synthGrating(varargin)
            % Constructor
            % See top comments for input documentation
            
            % adding path to emesim
            addpath( [ '..' filesep 'eme' ] );

            % inputs and defaults
            inputs = {  'discretization',   'none', ...
                        'units',            'um',   ...
                        'lambda',           'none', ...
                        'background_index', 1.0,    ...
                        'domain_size',      'none', ...
                        'period_vec',       'none', ...
                        'offset_vec',       'none', ...
                        'fill_top_vec',     'none', ...
                        'fill_bot_vec',     'none', ...
                        'optimal_angle',    'none', ...
                        'waveguide_index',  'none', ...
                        'waveguide_thicks', 'none', ...
                        'coupling_direction', 'down', ...
                        'data_directory',   'none', ...
                        'data_filename',    '', ...
                        'data_notes',       '', ...
                        'data_mode',        'new', ...
                        'num_par_workers',  'none', ...
                        'h_makeGratingCell', [] ...
                     }; 
%                         'ratio_vec',        'none', ...
%                         'fill_vec',         'none', ...
            obj.inputs = inputs;
            
            % first check whether to run code from fresh data or to load
            % previous results
            load_prev_result = false;   % defaults to starting fresh
            for ii = 1:2:length(varargin)
                if strcmp( varargin{ii}, 'data_mode' )
                    if strcmp( varargin{ii+1}, 'load' )
                        load_prev_result = true;
                    end
                end 
            end
            
            if ~load_prev_result
                % Starting a synth grating object from scratch
                fprintf('Starting a synth grating object from scratch\n\n');

                % parse inputs
                p = f_parse_varargin( inputs, varargin{:} );

                % save starting time
                obj.start_time = datestr( datetime('now'), 'yyyy_mm_dd HH_MM_SS ' );
                
                % set properties
%                 obj.fill_vec      = p.fill_vec;
%                 obj.ratio_vec     = p.ratio_vec;
                obj.fill_top_vec = p.fill_top_vec;
                obj.fill_bot_vec = p.fill_bot_vec;
                obj.period_vec    = p.period_vec;
                obj.offset_vec    = p.offset_vec;

                % set units
                obj.units.name  = p.units;
                switch( obj.units.name )
                    case 'm'
                        obj.units.scale = 1;
                    case 'mm'
                        obj.units.scale = 1e-3;
                    case 'um'
                        obj.units.scale = 1e-6;
                    case 'nm'
                        obj.units.scale = 1e-9;
                end

                obj.discretization      = p.discretization;
                obj.lambda              = p.lambda;
                obj.background_index    = p.background_index;
                obj.domain_size         = p.domain_size;
                obj.optimal_angle       = p.optimal_angle;

                obj.waveguide_index     = p.waveguide_index;
                obj.waveguide_thicks    = p.waveguide_thicks;

                if strcmp( p.coupling_direction, 'up') || strcmp( p.coupling_direction, 'down') 
                    % set coupling direction
                    obj.coupling_direction = p.coupling_direction;
                else
                    error('Error: input ''coupling_direction'' is not valid. Valid entries are ''up'' or ''down''. You entered ''%s''', p.coupling_direction);
                end

                obj.data_directory  = p.data_directory;
                obj.data_filename   = p.data_filename;
                obj.data_notes      = p.data_notes;
                obj.data_mode       = p.data_mode;
                
                obj.num_par_workers = p.num_par_workers;
                
                % default modesolver options (currently hardcoded)
                num_modes   = 20;
                BC          = 0;                    % 0 for PEC, 1 for PMC
                pml_options = [ 1, 200, 500, 2 ];   % [ yes/no, length in nm, strength, pml poly order ]
                obj.modesolver_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options );
                
                % set handle to grating cell making function
                if ~isempty( p.h_makeGratingCell )
                    obj.h_makeGratingCell = p.h_makeGratingCell;
                else
                    % default to the test cell at bottom of this class
                    obj.h_makeGratingCell = @makeGratingCell;
                end
                
                
            else
                % load previously run synth grating object
                fprintf('Loading a previously run synth grating object\n\n');
                
                % grab data directory and filename
                for ii = 1:2:length(varargin)
                    if strcmp( varargin{ii}, 'data_directory' )
                        data_directory = varargin{ii+1};
                    elseif strcmp( varargin{ii}, 'data_filename' )
                        data_filename = varargin{ii+1};
                    end 
                end
                
                obj             = obj.loadPreviousSweep(data_directory, data_filename);
                obj.data_mode   = 'load';
                
            end     % end if load previous result

        end     % end constructor()
       
        
        function obj = runParameterSweep( obj )
            % Runs full parameter sweep and saves all the data
            %
            % inputs:
            %   h_makeGratingCell
            %       type: Function handle
            %       desc: Handle to function that will instantiate and
            %             return a grating cell object, which can then be
            %             simulated in this parameter sweep code.
            %             MAKE SURE THE FUNCTION FOLLOWS THE EXACT CRITERIA
            %             OUTLINED IN THE TOP CLASS DOCUMENTATION
            %
            % example:
            %   obj = obj.runParameterSweep( @makeMyGratingCell )
            
            fprintf('Running parameter sweep...\n\n');
            
            % function handle to grating cell making function
            h_makeGratingCell = @obj.h_makeGratingCell;
         
            % extract some variables from the object
%             fill_vec    = obj.fill_vec;
%             ratio_vec   = obj.ratio_vec;
            fill_top_vec    = obj.fill_top_vec;
            fill_bot_vec    = obj.fill_bot_vec;
            period_vec      = obj.period_vec;
            offset_vec      = obj.offset_vec;

            % setup 4D tensors to save variable info
            % tensors have dimensions ( fill, ratio, period, offset )
%             [fill_tensor, ratio_tensor, period_tensor, offset_tensor] = ndgrid(fill_vec, ratio_vec, period_vec, offset_vec);
            [fill_top_tensor, fill_bot_tensor, period_tensor, offset_tensor] = ndgrid(fill_top_vec, fill_bot_vec, period_vec, offset_vec);
%             tensor_size         = size(fill_tensor);
            tensor_size         = size(fill_top_tensor);
            scatter_strengths   = zeros( tensor_size );
            directivities       = scatter_strengths;
            angles              = scatter_strengths;
            power_in            = scatter_strengths;
            power_rad_up        = scatter_strengths;
            power_rad_down      = scatter_strengths;
            
            
            % unwrap the tensors to make for easier looping, and thus
            % easier parallelization
%             fill_tensor         = fill_tensor(:);
%             ratio_tensor        = ratio_tensor(:);
            fill_top_tensor     = fill_top_tensor(:);
            fill_bot_tensor     = fill_bot_tensor(:);
            offset_tensor       = offset_tensor(:);
            period_tensor       = period_tensor(:);
            scatter_strengths   = scatter_strengths(:);
            directivities       = directivities(:);
            angles              = angles(:);
            power_in            = power_in(:);
            power_rad_up        = power_rad_up(:);
            power_rad_down      = power_rad_down(:);
            
            % run loops
%             num_loops   = length(fill_vec)*length(ratio_vec)*length(period_vec)*length(offset_vec);
            num_loops   = length(fill_top_vec)*length(fill_bot_vec)*length(period_vec)*length(offset_vec);
            
            
            % grab modesolver options
            num_modes   = obj.modesolver_opts.num_modes;
            BC          = obj.modesolver_opts.BC;        
            pml_options = obj.modesolver_opts.pml_options;
            
            % convert the object into a struct for the loop to use
            obj_copy = convertObjToStruct(obj);
            
            % start clock
            tic;
            
            
            % init parallel pool
            
            % Taken from the BU SCC documentation:
            % Especially important for running multiple batch jobs
            % Without this procedure, some batch jobs may fail
            % redirects ~/.matlab PCT temp files to system's TMPDIR on compute
            % node to avoid inter-node (compute node <--> login node) I/O
            myCluster = parcluster('local');                        % cores on compute node are "local"
            if getenv('ENVIRONMENT')                                % true if this is a batch job
                myCluster.JobStorageLocation = getenv('TMPDIR');    % points to TMPDIR
            end

            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if ~isempty(poolobj)
                % shut down previously made parallel pool
                delete(gcp('nocreate'));
            end
            parpool(myCluster, obj.num_par_workers);
            
            
            
            parfor ii = 1:num_loops
                
                fprintf('Running loop %i of %i\n', ii, num_loops);
                
                % grab some parameters
                period          = period_tensor(ii);
%                 fill            = fill_tensor(ii);
%                 ratio           = ratio_tensor(ii);
                fill_top        = fill_top_tensor(ii);
                fill_bot        = fill_bot_tensor(ii);
                offset_ratio    = offset_tensor(ii);
                
                % make grating cell
%                 Q = h_makeGratingCell( obj_copy, period, fill, ratio, offset_ratio );
                Q = h_makeGratingCell( obj_copy, period, fill_top, fill_bot, offset_ratio );
                
                % run simulation
                Q = Q.runSimulation( num_modes, BC, pml_options );
                
                % save parameters
                if strcmp(obj.coupling_direction, 'up')
                    % coupling direction is up
                    directivities(ii)       = Q.directivity;
                    scatter_strengths(ii)   = Q.alpha_up;
                    angles(ii)              = Q.max_angle_up;
                else
                    % coupling direction is down
                    directivities(ii)       = 1/Q.directivity;
                    scatter_strengths(ii)   = Q.alpha_down;
                    angles(ii)              = Q.max_angle_down;
                end
                
                power_in(ii)        = Q.P_in;
                power_rad_up(ii)    = Q.P_rad_up;
                power_rad_down(ii)  = Q.P_rad_down;
                
            end     % end for ii = 1:num_loops
            toc;
            
            % free the parallel pool
            delete(gcp('nocreate'));
            
            % reshape the unwrapped tensors back into tensor form
%             fill_tensor         = reshape( fill_tensor, tensor_size );
%             ratio_tensor        = reshape( ratio_tensor, tensor_size );
            offset_tensor       = reshape( offset_tensor, tensor_size );
            period_tensor       = reshape( period_tensor, tensor_size );
            fill_top_tensor     = reshape( fill_top_tensor, tensor_size );
            fill_bot_tensor     = reshape( fill_bot_tensor, tensor_size );
            scatter_strengths   = reshape( scatter_strengths, tensor_size );
            directivities       = reshape( directivities, tensor_size );
            angles              = reshape( angles, tensor_size );
            power_in            = reshape( power_in, tensor_size );
            power_rad_up        = reshape( power_rad_up, tensor_size );
            power_rad_down      = reshape( power_rad_down, tensor_size );
            
            % save all data to a mat file
%             sweep_results = struct( ...  % 'fill_tensor', fill_tensor, ...
% %                                     'ratio_tensor', ratio_tensor, ...
%                                     'fill_top_tensor', fill_top_tensor, ...
%                                     'fill_bot_tensor', fill_bot_tensor, ...
%                                     'offset_tensor', offset_tensor, ...
%                                     'period_tensor', period_tensor, ...
%                                     'scatter_strengths', scatter_strengths, ...
%                                     'directivities', directivities, ...
%                                     'angles', angles, ...
%                                     'power_in', power_in, ...
%                                     'power_rad_up', power_rad_up, ...
%                                     'power_rad_down', power_rad_down );                     
            sweep_results = struct( 'fill_top_tensor', fill_top_tensor, ...
                                    'fill_bot_tensor', fill_bot_tensor, ...
                                    'offset_tensor', offset_tensor, ...
                                    'period_tensor', period_tensor, ...
                                    'scatter_strengths', scatter_strengths, ...
                                    'directivities', directivities, ...
                                    'angles', angles, ...
                                    'power_in', power_in, ...
                                    'power_rad_up', power_rad_up, ...
                                    'power_rad_down', power_rad_down );
                                
            % store sweep results to synthgrating object
            obj.sweep_results = sweep_results;
            
            
            full_filename   = [ obj.data_directory, filesep, obj.start_time, obj.data_filename, '.mat' ];
            
            fprintf('Saving data to directory %s, filename ''%s''...\n', obj.data_directory, full_filename);
%             save(full_filename,'synth_obj');
            obj.saveToStruct( full_filename );
            fprintf('...done\n\n');
             
%             % DEBUG testing waitbar
%             for ii = 1:num_loops
%                waitbar(ii/num_loops, h_waitbar, sprintf('Loop %i of %i', ii, num_loops));
%                pause(1);
%             end
            
%             % close the waiting bar
%             delete( h_waitbar );
            
            fprintf('\n...done running parameter sweep\n\n');
            
        end         

        
        function obj = loadPreviousSweep(obj, data_directory, data_filename)
            % This function is for loading a previous synthGrating object
            % that already has performed a parameter sweep.
            
            % this loads the "sweep_obj" struct
            load( [ data_directory, filesep, data_filename ] );
            
            fields = fieldnames( sweep_obj );
            
            for ii = 1:length(fields)
                % overwrite the current object's data (if that field
                % exists)
                if isprop( obj, fields{ii} ) 
                    obj.(fields{ii}) = sweep_obj.(fields{ii});
                end
            end
            
            % TEMPORARY, write handle to grating drawing function if it
            % doesn't exist
            if isempty( obj.h_makeGratingCell )
                obj.h_makeGratingCell = @makeGratingCell;
            end

        end
        
        
        function saveToStruct(obj, filename)
            % Saves all current properties of this object to a structure,
            % and then to a .mat file
            
            sweep_obj = obj.convertObjToStruct();
            save(filename, 'sweep_obj');

        end
        
        
        function obj_as_struct = convertObjToStruct(obj)
            % converts the current object to a struct that holds the
            % object's properties
            
            props = properties(obj);
     
            obj_as_struct = struct();
            for p = 1:numel(props)
                if strcmp( props{p}, 'h_makeGratingCell' )
                    % convert function handle to string
                    obj_as_struct.(props{p}) = func2str( obj.(props{p}) );
                else
                    obj_as_struct.(props{p}) = obj.(props{p});
                end
            end
            
        end
            

        function [obj, u] = fiberModeGaussian(obj, w0, zvec, xvec, theta, d0, nclad)
        % somewhat adapted from Cale's code
        %
        % Generate Gaussian-beam mode profile at a plane through y = d0 at angle of
        % theta (in degrees)
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
        %   zvec
        %       type: double, array
        %       desc: out of plane coordinates... currently not used
        %
        %   xvec
        %       type: double, array
        %       desc: coordinates along direction of propagation of
        %             grating, in units 'units'
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
        %
        %
        % Outputs: 
        %   u
        %       type: double, array
        %       desc: returned slice of gaussian beam, normalized to total
        %             power


            % Constants, in units of meters
            lambda  = obj.lambda * obj.units.scale / nclad;     % wavelength in cladding, units m
            k0      = 2*pi/lambda;                              % 1/m
            w0      = w0 * obj.units.scale;                     % [meters] radius
            d0      = d0 * obj.units.scale;                     % [meters] offset
            
            % Convert to radians
            theta = (pi/180)*theta;
            
            % Scale coordinates
            xvec = xvec * obj.units.scale;                                              % units m
            yvec = xvec;
            zvec = zvec * obj.units.scale;                                              % units m
            
            % try just plotting this slice of data
            xprime = xvec.*cos(-theta) + d0*sin(-theta);
            zprime = -xvec.*sin(-theta) + d0*cos(-theta);

            % b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)
            b = k0*w0^2/2;                                                                                   

            % Equation (5.2) in Haus [1/meters]
            u00_slice =   1j .* sqrt(k0*b/pi) .* ( 1./(zprime + 1j.*b) ).*...
                exp( -1j.*k0.*( xprime.^2 )./( 2*(zprime + 1j.*b) ) );     
            
            % normalize the slice to intensity
            dx          = obj.discretization * obj.units.scale;                 % disc. in m
            u00_slice   = u00_slice/sqrt( dx * sum( abs( u00_slice ).^2 ) );
            
            % return and save data
            u       = u00_slice;
            obj.u   = u;

        end     % end fiberModeGaussian()

        function obj = synthesizeUniformGrating(obj, angle, MFD)
            % THIS FUNCTION IS DEPRECATED
            %
            % "Synthesizes" a uniform grating at the desired angle 
            % more of a simple test function to see whether the two level
            % grating simulation, parameter sweep, and gaussian overlap are
            % working.
            %
            % inputs:
            %   angle in deg
            %   MFD in units 'units'
            
            % number of periods
            N_periods = 20;
            
%             % pick the parameters that resulted in output angle closest to
%             % desired angle
%             [ ~, i_min ]                                = min( abs(obj.sweep_results.angles(:) - angle) );
%             [ i_fill, i_ratio, i_period, i_offset ]     = ind2sub( size(obj.sweep_results.angles), i_min );
%             
%             % DEBUG show results
%             % tensors have dimensions ( fill, ratio, period, offset )
%             fill    = obj.fill_vec( i_fill );
%             ratio   = obj.ratio_vec( i_ratio );
%             offset  = obj.offset_vec( i_offset );
%             period  = obj.period_vec( i_period );
%             angle               = obj.sweep_results.angles( i_fill, i_ratio, i_period, i_offset );
%             directivity         = obj.sweep_results.directivities( i_fill, i_ratio, i_period, i_offset );
%             scatter_strength    = obj.sweep_results.scatter_strengths( i_fill, i_ratio, i_period, i_offset );
            
            % new version, pick parameters that result in output angle
            % being within tolerance of a few deg
            % first unwrap all the variables
            fills               = obj.sweep_results.fill_tensor(:);
            ratios              = obj.sweep_results.ratio_tensor(:);
            offsets             = obj.sweep_results.offset_tensor(:);
            periods             = obj.sweep_results.period_tensor(:); 
            angles              = obj.sweep_results.angles(:);
            directivities       = obj.sweep_results.directivities(:);
            scatter_strengths   = obj.sweep_results.scatter_strengths(:);
            
            % now pick the designs with the angles within tolerance
            tol                 = 5;                                                % angle tolerance, in degrees
            indx_best_angles    = abs(angles - angle) <= tol;                       % indexes of grating cells that fit in this tolerance
            fills               = fills( indx_best_angles );
            ratios              = ratios( indx_best_angles );
            offsets             = offsets( indx_best_angles );
            periods             = periods( indx_best_angles );
            angles              = angles( indx_best_angles );
            directivities       = directivities( indx_best_angles );
            scatter_strengths      = scatter_strengths( indx_best_angles );
            
% %             % now pick the designs with the scattering strength within tolerance
% %             % scattering strength tolerance, in percent from maximum
% %             tol                 = 5;
% %             % indexes of grating cells that fit in this tolerance
% %             indx_best_scatter   = abs(scatter_strengths - max(scatter_strengths(:)))./max(scatter_strengths(:)) <= tol/100;                       
% %             fills               = fills( indx_best_scatter );
% %             ratios              = ratios( indx_best_scatter );
% %             offsets             = offsets( indx_best_scatter );
% %             periods             = periods( indx_best_scatter );
% %             angles              = angles( indx_best_scatter );
% %             directivities       = directivities( indx_best_scatter );
% %             scatter_strengths      = scatter_strengths( indx_best_scatter );
%             
            % now pick the design with the best directivity
            [max_directivity, indx_best_dir] = max( directivities );
            fill            = fills( indx_best_dir );
            ratio           = ratios( indx_best_dir );
            offset          = offsets( indx_best_dir );
            period          = periods( indx_best_dir );
            angle_actual    = angles( indx_best_dir );
            scatter_strength   = scatter_strengths( indx_best_dir );
            
%             max_directivity = max(obj.sweep_results.directivities(:));

            % use these parameters to build a full grating
            % gonna use an emeSim to build this grating
            % Set Up Simulation
            % note that emeSim uses 'z' as propagation direction and 'x'
            % as transverse (synthGrating uses 'x' and 'y' respectively)
            % and units are in um
            um          = 1e6;
            dx          = obj.discretization * obj.units.scale * um;                % in um
            dz          = 5e-3;                                                     % in um
            pol         = 0;                                                        % 0 for TE, 1 for TM
            z_in        = 1.5;                                                      % length of input section of waveguide
            xf          = obj.domain_size(1) * obj.units.scale * um;                % in um
            zf          = N_periods*period*obj.units.scale*um + z_in;               % in um
            lambda_um   = obj.lambda * obj.units.scale * um;                        % wl in um
            eme_obj = emeSim(   'discretization', [dx dz], ...
                                'pml', 0.2, ...
                                'domain', [xf zf], ...
                                'backgroundIndex', obj.background_index, ...
                                'wavelengthSpectrum', [lambda_um lambda_um 0.1], ...
                                'debug', 'no',...                   
                                'polarization', pol );
                            
            % TRICK - i can use the twoLeveLgratingcell to build the
            % dielectric for the emeSim
            % make grating cell, in units of um
            domain      = [ xf, period * obj.units.scale * um ];
            gratingcell = c_twoLevelGratingCell(  'discretization', [ dx, dz ], ...
                                        'units', 'um', ...
                                        'lambda', lambda_um, ...
                                        'domain_size', domain, ...
                                        'background_index', obj.background_index );

            % draw cell
            % draw two levels using two level builder function
            wg_thick_um     = obj.waveguide_thicks * obj.units.scale * um;     % in um
            wg_min_y        = [ domain(1)/2, domain(1)/2-wg_thick_um(1) ];
            wgs_duty_cycles = [ fill, fill*ratio ];
            wgs_offsets     = [ 0, offset*period*obj.units.scale*um ];
            gratingcell     = gratingcell.twoLevelBuilder(  wg_min_y, wg_thick_um, obj.waveguide_index, ...
                                                            wgs_duty_cycles, wgs_offsets );
                                                        
                                                        
            % create cell dielectric
            cell_diel = gratingcell.N;
            
            
            % stitch together all the cells
            diel            = repmat( cell_diel, 1, N_periods );
            eme_obj.diel    = diel;
            
            
            % draw the input waveguide section and stitch that as well
            domain          = [ xf, z_in ];
            gratingcell_in  = c_twoLevelGratingCell(  'discretization', [ dx, dz ], ...
                                        'units', 'um', ...
                                        'lambda', lambda_um, ...
                                        'domain_size', domain, ...
                                        'background_index', obj.background_index );
            % draw input wg
            wg_thick_um     = obj.waveguide_thicks * obj.units.scale * um;     % in um
            wg_min_y        = [ domain(1)/2, domain(1)/2-wg_thick_um(1) ];
            gratingcell_in  = gratingcell_in.twoLevelBuilder(  wg_min_y, wg_thick_um, obj.waveguide_index, ...
                                                            [ 1, 1 ], [ 0, 0 ] );
            eme_obj.diel    = [ gratingcell_in.N, eme_obj.diel ];
            
%             % DEBUG plot the diel
%             eme_obj.plotDiel();
            
            % run EME sim
            % Converts the dielectric distribution into layers for eigen mode expansion
            eme_obj = eme_obj.convertDiel();   
            % Runs simulation
            eme_obj = eme_obj.runSimulation('plotSource','yes');      
            % compute fiber overlap
            eme_obj = eme_obj.fiberOverlap( 'zOffset', 0:.1:12,...
                                            'angleVec', -45:1:45,...
                                            'MFD', MFD * obj.units.scale * um,...
                                            'overlapDir', obj.coupling_direction);
                                        
            % DEBUG show results
            gratingUI(eme_obj);

        end
        
        function obj = synthesizeGaussianGrating(obj, angle, MFD)
            % Synthesizes a grating that is mode-matched to fiber gaussian
            % mode
            %
            % Inputs:
            %   angle
            %       type: double, scalar
            %       desc: angle of fiber from normal
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter
            
            % number of cells
            n_cells = 20;
            
            % generate x coordinates for the gaussian mode
            % must be large enough to fit all cells + mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiberModeGaussian(    w0, zvec, xvec,...
                                                    angle, d0, obj.background_index );
            
            % calculate desired scattering strength vs. x
            integral_u  = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
            alpha_des   = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );             % in units 1/m
            alpha_des   = alpha_des * obj.units.scale;                                  % in units 1/units

%             % DEBUG plot u
%             figure;
%             plot( xvec, abs(u) );
%             xlabel(['x (' obj.units.name ')']); title('DEBUG slice of gaussian');
%             makeFigureNice();
% 
%             % DEBUG plot integral_u
%             figure;
%             plot( xvec, integral_u );
%             xlabel(['x (' obj.units.name ')']); title('DEBUG integral of gaussian^2');
%             makeFigureNice();
%             
            % DEBUG plot alpha desired
            figure;
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice();

            % DEBUG unwrap and see what variable ranges have been simulated
            % first unwrap all the variables
%             fills               = obj.sweep_results.fill_tensor(:);
%             ratios              = obj.sweep_results.ratio_tensor(:);
            fill_tops           = obj.sweep_results.fill_top_tensor(:);
            fill_bots           = obj.sweep_results.fill_bot_tensor(:);
            offsets             = obj.sweep_results.offset_tensor(:);
            periods             = obj.sweep_results.period_tensor(:); 
            angles              = obj.sweep_results.angles(:);
            directivities       = obj.sweep_results.directivities(:);
            scatter_strengths   = obj.sweep_results.scatter_strengths(:);
            
            % DEBUG sort and plot all the simulated angles
            angles_sorted = sort( angles );
            figure; 
            plot( 1:length(angles), angles_sorted );
            makeFigureNice();
            title('DEBUG all simulated angles, sorted in ascending order');
            
            % DEBUG sort and plot all the simulated directivities
            directivities_sorted = sort( directivities );
            figure; 
            plot( 1:length(directivities), directivities_sorted );
            makeFigureNice();
            title('DEBUG all simulated directivities, sorted in ascending order');
            
            % DEBUG sort and plot all the simulated scatter_strengths
            scatter_strengths_sorted = sort( scatter_strengths );
            figure; 
            plot( 1:length(directivities), scatter_strengths_sorted );
            makeFigureNice();
            title('DEBUG all simulated scatter strengths, sorted in ascending order');

            % init 2D data variables
            % dimensions fill x ratio
            chosen_angles           = squeeze( zeros( size( obj.sweep_results.angles( :, :, 1, 1 ) ) ) );
            chosen_periods          = zeros( size(chosen_angles) );
            chosen_directivities    = zeros( size(chosen_angles) );
            chosen_scatter_str      = zeros( size(chosen_angles) ); 
            chosen_offsets          = zeros( size(chosen_angles) );
            chosen_closest_angles   = zeros( size(chosen_angles) );     % DEBUG, what the angle would be without changing offset
            
            % Synthesis loop
            % OLD VERSION, using fill and ratio
%             for i_fill = 1:length(obj.fill_vec)
%                 % for each fill
%                 for i_ratio = 1:length(obj.ratio_vec)
%                     % for each ratio
%                     
%                     % pick period with angle closest to desired
%                     % tensors have dimensions ( fill, ratio, period, offset )
%                     angles_per_fill_ratio   = squeeze( obj.sweep_results.angles( i_fill, i_ratio, :, : ) );         % dimensiosn period x offset
%                     dirs_per_fill_ratio     = squeeze( obj.sweep_results.directivities( i_fill, i_ratio, :, : ) );  % dimensiosn period x offset
% 
% %                     % index of angle and period closest to desired
%                     [ ~, angle_indx ]                       = min( abs(angles_per_fill_ratio(:) - angle) );
%                     [ i_period, i_offset_closest_angle ]    = ind2sub( size(angles_per_fill_ratio), angle_indx );
% 
% 
%                     % instead let's try a multi objective optimization
%                     % minimizing the angle and maximizing the directivity
%                     % together
%                     min_angle_merit         = abs( angles_per_fill_ratio - angle )/abs(angle);
%                     max_dir_merit           = 1 - 10*log10(dirs_per_fill_ratio) / max( 10*log10(dirs_per_fill_ratio(:)) );
%                     merit                   = 10*min_angle_merit + max_dir_merit;      % total merit function to minimize
%                     [ ~, indx_best_merit ]  = min( merit(:) );
%                     [ i_period, i_offset ]  = ind2sub( size(merit), indx_best_merit );
%                     
%                     % save the chosen variables
% %                     chosen_angles( i_fill, i_ratio )        = angles_per_fill_ratio( i_period, i_offset );
%                     chosen_angles( i_fill, i_ratio )        = obj.sweep_results.angles( i_fill, i_ratio, i_period, i_offset );
% %                     chosen_directivities( i_fill, i_ratio ) = max_dir;
%                     chosen_directivities( i_fill, i_ratio ) = obj.sweep_results.directivities( i_fill, i_ratio, i_period, i_offset );
%                     chosen_periods( i_fill, i_ratio )       = obj.sweep_results.period_tensor( i_fill, i_ratio, i_period, i_offset );
%                     chosen_offsets( i_fill, i_ratio )       = obj.sweep_results.offset_tensor( i_fill, i_ratio, i_period, i_offset );
%                     chosen_scatter_str( i_fill, i_ratio )   = obj.sweep_results.scatter_strengths( i_fill, i_ratio, i_period, i_offset );
%                     chosen_closest_angles( i_fill, i_ratio) = obj.sweep_results.angles( i_fill, i_ratio, i_period, i_offset_closest_angle );
%                     
%                     
%                 end
%             end
            
            % Synthesis loop, using fill top and fill bot
            for i_fill_top = 1:length(obj.fill_top_vec)
                % for each fill
                for i_fill_bot = 1:length(obj.fill_bot_vec)
                    % for each ratio
                    
                    % pick period with angle closest to desired
                    % tensors have dimensions ( fill, ratio, period, offset )
                    angles_per_fill_ratio   = squeeze( obj.sweep_results.angles( i_fill_top, i_fill_bot, :, : ) );         % dimensiosn period x offset
                    dirs_per_fill_ratio     = squeeze( obj.sweep_results.directivities( i_fill_top, i_fill_bot, :, : ) );  % dimensiosn period x offset

%                     % index of angle and period closest to desired
                    [ ~, angle_indx ]                       = min( abs(angles_per_fill_ratio(:) - angle) );
                    [ i_period, i_offset_closest_angle ]    = ind2sub( size(angles_per_fill_ratio), angle_indx );


                    % instead let's try a multi objective optimization
                    % minimizing the angle and maximizing the directivity
                    % together
                    min_angle_merit         = abs( angles_per_fill_ratio - angle )/abs(angle);
                    max_dir_merit           = 1 - 10*log10(dirs_per_fill_ratio) / max( 10*log10(dirs_per_fill_ratio(:)) );
                    merit                   = 10*min_angle_merit + max_dir_merit;      % total merit function to minimize
                    [ ~, indx_best_merit ]  = min( merit(:) );
                    [ i_period, i_offset ]  = ind2sub( size(merit), indx_best_merit );
                    
                    % save the chosen variables
%                     chosen_angles( i_fill, i_ratio )        = angles_per_fill_ratio( i_period, i_offset );
                    chosen_angles( i_fill_top, i_fill_bot )        = obj.sweep_results.angles( i_fill_top, i_fill_bot, i_period, i_offset );
%                     chosen_directivities( i_fill, i_ratio ) = max_dir;
                    chosen_directivities( i_fill_top, i_fill_bot ) = obj.sweep_results.directivities( i_fill_top, i_fill_bot, i_period, i_offset );
                    chosen_periods( i_fill_top, i_fill_bot )       = obj.sweep_results.period_tensor( i_fill_top, i_fill_bot, i_period, i_offset );
                    chosen_offsets( i_fill_top, i_fill_bot )       = obj.sweep_results.offset_tensor( i_fill_top, i_fill_bot, i_period, i_offset );
                    chosen_scatter_str( i_fill_top, i_fill_bot )   = obj.sweep_results.scatter_strengths( i_fill_top, i_fill_bot, i_period, i_offset );
                    chosen_closest_angles( i_fill_top, i_fill_bot) = obj.sweep_results.angles( i_fill_top, i_fill_bot, i_period, i_offset_closest_angle );
                    
                    
                end
            end
            
            % IF YOU WANT TO PLOT THESE AGAIN, replace obj.ratio_vec with
            % obj.fill_bot_vec, and obj.fill_vec with obj.fill_top_vec
%             % DEBUG plot the 2D design spaces
%             % chosen angles
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, chosen_angles );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for angles vs. fill and ratio');
%             colorbar;
%             % chosen directivities
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, 10*log10(chosen_directivities) );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for directivities (dB) vs. fill and ratio');
%             colorbar;
%             % chosen scatter strengths
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, chosen_scatter_str );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for scatter strengths vs. fill and ratio');
%             colorbar;
%             % chosen periods
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, chosen_periods );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for periods vs. fill and ratio');
%             colorbar;
%             % chosen offsets
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, chosen_offsets );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for offsets vs. fill and ratio');
%             colorbar;
%             % DEBUG chosen closest angles 
%             figure;
%             imagesc( obj.ratio_vec, obj.fill_vec, chosen_closest_angles );
%             xlabel('ratios'); ylabel('fill');
%             set(gca, 'ydir', 'normal');
%             title('DEBUG plot of chosen 2D design space for CLOSEST angles vs. fill and ratio');
%             colorbar;

%             % Synthesis curve v1
%             % For each fill value, pick the design with the highest
%             % directivity.
%             % plot the other resulting parameters
%             
%             % init the data saving variables
%             best_dir_v_fill     = zeros( size(obj.fill_vec) );
%             best_period_v_fill  = best_dir_v_fill;
%             best_angle_v_fill   = best_dir_v_fill;
%             best_offset_v_fill  = best_dir_v_fill;
%             best_scatter_v_fill = best_dir_v_fill;
%             best_ratio_v_fill   = best_dir_v_fill;
%             
%             for ii = 1:length( obj.fill_vec )
%                 
%                 % grab the directivities for this fill value
%                 cur_dir = chosen_directivities( ii, : );
%                 
%                 % find max directivity
%                 [ max_dir, indx_max_dir ] = max( cur_dir(:) );
%                 
%                 % save stuff
%                 best_dir_v_fill(ii)     = max_dir;
%                 best_period_v_fill(ii)  = chosen_periods( ii, indx_max_dir );
%                 best_angle_v_fill(ii)   = chosen_angles( ii, indx_max_dir );
%                 best_offset_v_fill(ii)  = chosen_offsets( ii, indx_max_dir );
%                 best_scatter_v_fill(ii) = chosen_scatter_str( ii, indx_max_dir );
%                 best_ratio_v_fill(ii)   = obj.ratio_vec( indx_max_dir );
%                 
%             end
%             
%             % plot these results
%             % directionality
%             figure;
%             plot( obj.fill_vec, 10*log10(best_dir_v_fill), '-o' );
%             xlabel('fill ratio'); ylabel('dB');
%             title('best directionality vs. fill ratio');
%             makeFigureNice();
%             % ratio
%             figure;
%             plot( obj.fill_vec, best_ratio_v_fill, '-o' );
%             xlabel('fill ratio'); ylabel('layer ratio');
%             title('Layer ratio for best directionality vs. fill ratio');
%             makeFigureNice();
            
            % Synthesis curve v2
            % For each ratio value, pick the design with the highest
            % directivity.
            % plot the other resulting parameters
            
            
            
            % -!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-
            % STUFF BELOW HERE NEES TO BE REDONE
            
            
            % init the data saving variables
            best_dir_v_ratio     = zeros( size(obj.ratio_vec) );
            best_period_v_ratio  = best_dir_v_ratio;
            best_angle_v_ratio   = best_dir_v_ratio;
            best_offset_v_ratio  = best_dir_v_ratio;
            best_scatter_v_ratio = best_dir_v_ratio;
            best_fill_v_ratio    = best_dir_v_ratio;
            
            for ii = 1:length( obj.ratio_vec )
                
                % grab the directivities for this ratio value
                cur_dir = chosen_directivities( :, ii );
                
                % find max directivity
                [ max_dir, indx_max_dir ] = max( cur_dir(:) );
                
                % save stuff
                best_dir_v_ratio(ii)     = max_dir;
                best_period_v_ratio(ii)  = chosen_periods( indx_max_dir, ii );
                best_angle_v_ratio(ii)   = chosen_angles( indx_max_dir, ii  );
                best_offset_v_ratio(ii)  = chosen_offsets( indx_max_dir, ii  );
                best_scatter_v_ratio(ii) = chosen_scatter_str( indx_max_dir, ii  );
                best_fill_v_ratio(ii)    = obj.fill_vec( indx_max_dir );
                
            end
            
            % plot these results
            % directionality
            figure;
            plot( obj.ratio_vec, 10*log10(best_dir_v_ratio), '-o' );
            xlabel('layer ratio'); ylabel('dB');
            title('best directionality vs. layer ratio');
            makeFigureNice();
            % fill
            figure;
            plot( obj.ratio_vec, best_fill_v_ratio, '-o' );
            xlabel('layer ratio'); ylabel('fill ratio');
            title('Fill ratio for best directionality vs. layer ratio');
            makeFigureNice();
            % angle
            figure;
            plot( obj.ratio_vec, best_angle_v_ratio, '-o' );
            xlabel('layer ratio'); ylabel('angle');
            title('Angle for best directionality vs. layer ratio');
            makeFigureNice();
            % period
            figure;
            plot( obj.ratio_vec, best_period_v_ratio, '-o' );
            xlabel('layer ratio'); ylabel( obj.units.name );
            title('Period for best directionality vs. layer ratio');
            makeFigureNice();
            % offset
            figure;
            plot( obj.ratio_vec, best_offset_v_ratio, '-o' );
            xlabel('layer ratio'); ylabel('offset ratio');
            title('Offset ratio for best directionality vs. layer ratio');
            makeFigureNice();
            
%             % DEBUG
%             % unfortunately I fucked up and didn't simulate the scattering
%             % strength correctly, so as a workaround I have to re-simulate
%             % these gratings again
% %             tic;
%             
%             % start parpool
%             poolobj = gcp('nocreate'); % If no pool, do not create new one.
%             if ~isempty(poolobj)
%                 % shut down previously made parallel pool
%                 delete(gcp('nocreate'));
%             end
%             parpool('local', 4);
% 
%             n_loops     = length(obj.ratio_vec);
%             parfor ii = 1:length( obj.ratio_vec )
%                
%                 fprintf('DEBUG running re-sim loop %i of %i\n', ii, n_loops );
%                 
%                 % simulate grating
%                 GC = obj.h_makeGratingCell( obj.convertObjToStruct(), best_period_v_ratio(ii), best_fill_v_ratio(ii), ...
%                                             obj.ratio_vec(ii), best_offset_v_ratio(ii) );
%                                         
%                 % run simulation
%                 GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
%                 
%                 % save scattering strength
%                 if strcmp(obj.coupling_direction, 'up')
%                     best_scatter_v_ratio(ii) = GC.alpha_up;
%                 else
%                     best_scatter_v_ratio(ii) = GC.alpha_down;
%                 end
%                 
% %                 toc;
%                 
%             end     % end for ii = 1:length(obj.ratio_vec)
            
            % plot scatter
            figure;
            plot( obj.ratio_vec, best_scatter_v_ratio, '-o' );
            xlabel('layer ratio'); ylabel('scatter strength');
            title('Scatter strength (\alpha) for best directionality vs. layer ratio');
            makeFigureNice();
            
            % LAST PART of the synthesis
            % now for each cell, pick parameters that give closest
            % scattering strength
            
            % DEBUG cut off the values beyond layer ratio of 0.8 lol
            indx_below_0d8 = obj.ratio_vec <= 0.8;
            % save stuff
            best_dir_v_ratio       = best_dir_v_ratio( indx_below_0d8 );
            best_period_v_ratio    = best_period_v_ratio( indx_below_0d8 );
            best_angle_v_ratio     = best_angle_v_ratio( indx_below_0d8 );
            best_offset_v_ratio    = best_offset_v_ratio( indx_below_0d8 );
            best_scatter_v_ratio   = best_scatter_v_ratio( indx_below_0d8 );
            best_fill_v_ratio      = best_fill_v_ratio( indx_below_0d8 );
            best_ratio_vec         = obj.ratio_vec( indx_below_0d8 );

            % first pick starting point for gaussian
            xstart          = -MFD/2;
            [~, indx_x]     = min( abs(xvec - xstart) );
            cur_x           = xvec( indx_x );
            
            % save data to these variables
            max_directivities_synth    = zeros( 1, n_cells );
            fills_synth                = zeros( 1, n_cells );
            ratios_synth               = zeros( 1, n_cells );
            offsets_synth              = zeros( 1, n_cells );
            periods_synth              = zeros( 1, n_cells );
            angles_synth               = zeros( 1, n_cells );
            scatter_strengths_synth    = zeros( 1, n_cells );
            
            % let's try normalizing the alphas, because I know that they
            % aren't on the same order of magnitude right now.
            alpha_des_norm              = alpha_des./max(alpha_des(:));
            best_scatter_v_ratio_norm   = best_scatter_v_ratio./max(best_scatter_v_ratio(:));
            best_scatter_v_ratio_norm   = best_scatter_v_ratio_norm - min(best_scatter_v_ratio_norm);   % DEBUG what if i remove the DC
            
            % DEBUG plot alpha des and best scatter
            figure;
            plot( 1:length(alpha_des_norm), alpha_des_norm );
            title('desired normalized \alpha'); 
            makeFigureNice();
            
            figure;
            plot( 1:length( best_scatter_v_ratio_norm ),best_scatter_v_ratio_norm, '-o');
            title('\alpha''s to choose from');
            makeFigureNice();
            
            
            for ii = 1:n_cells
                
%                 fprintf('%i\n', ii); % DEBUG
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter             = alpha_des_norm(indx_x);                                    % desired alpha
                [~, indx_best_scatter]  = min( abs(best_scatter_v_ratio_norm - des_scatter) );      % index of closest scatter design         
                
                % save parameters
                max_directivities_synth(ii) = best_dir_v_ratio( indx_best_scatter );
                fills_synth(ii)             = best_fill_v_ratio( indx_best_scatter );
                ratios_synth(ii)            = best_ratio_vec( indx_best_scatter );
                offsets_synth(ii)           = best_offset_v_ratio( indx_best_scatter );
                periods_synth(ii)           = best_period_v_ratio( indx_best_scatter );
                angles_synth(ii)            = best_angle_v_ratio( indx_best_scatter );
                scatter_strengths_synth(ii) = best_scatter_v_ratio( indx_best_scatter );
                
                % move onto next
                cur_x       = cur_x + periods_synth(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                
            end     % end for ii = 1:ncells
            
            % plot the synthesized design
            % scatter strength
            figure;
            plot( 1:n_cells, scatter_strengths_synth, '-o' );
            xlabel('cell #'); ylabel(['\alpha (1/' obj.units.name ')']);
            title('Synthesized scattering strengths v cell');
            makeFigureNice();
            % period
            figure;
            plot( 1:n_cells, periods_synth, '-o' );
            xlabel('cell #'); ylabel( obj.units.name );
            title('Synthesized periods v cell');
            makeFigureNice();
            % layer ratio
            figure;
            plot( 1:n_cells, ratios_synth, '-o' );
            xlabel('cell #'); ylabel( 'layer ratio' );
            title('Synthesized layer ratio v cell');
            makeFigureNice();
            % fill ratio
            figure;
            plot( 1:n_cells, fills_synth, '-o' );
            xlabel('cell #'); ylabel( 'fill ratio' );
            title('Synthesized fill ratio v cell');
            makeFigureNice();
            % offset ratio
            figure;
            plot( 1:n_cells, offsets_synth, '-o' );
            xlabel('cell #'); ylabel( 'offset ratio' );
            title('Synthesized offset ratio v cell');
            makeFigureNice();
            % angle
            figure;
            plot( 1:n_cells, angles_synth, '-o' );
            xlabel('cell #'); ylabel( 'deg' );
            title('Synthesized angle v cell');
            makeFigureNice();
            
            % NOW to verify the design
            % Run it in EME
            % Set Up Simulation
            % note that emeSim uses 'z' as propagation direction and 'x'
            % as transverse (synthGrating uses 'x' and 'y' respectively)
            % and units are in um
            um          = 1e6;
            dx          = obj.discretization * obj.units.scale * um;                % in um
            dz          = 5e-3;                                                     % in um
            pol         = 0;                                                        % 0 for TE, 1 for TM
            z_in        = 1.5;                                                      % length of input section of waveguide
            xf          = obj.domain_size(1) * obj.units.scale * um;                % in um
            zf          = sum( periods_synth )*obj.units.scale*um + z_in;           % in um
            lambda_um   = obj.lambda * obj.units.scale * um;                        % wl in um
            eme_obj     = emeSim(   'discretization', [dx dz], ...
                                'pml', 0.2, ...
                                'domain', [xf zf], ...
                                'backgroundIndex', obj.background_index, ...
                                'wavelengthSpectrum', [lambda_um lambda_um 0.1], ...
                                'debug', 'no',...                   
                                'polarization', pol );
            diel        = eme_obj.diel;
            % grab emeSim coordinates
            z_coords_eme    = eme_obj.domain.z;
            cur_z           = z_coords_eme(1);          % current z coordinate
            
            % draw the input waveguide section
            % using the trick that i can write and return the index from
            % the two level grating cell
            % first override the discretization
            obj_as_struct                   = obj.convertObjToStruct();
            obj_as_struct.discretization    = [ dx, dz ] / ( um * obj.units.scale );
            % now make the grating cell
            gratingcell_in  = obj.h_makeGratingCell( obj_as_struct, ...
                                                     z_in/(um*obj.units.scale), ...
                                                     1.0, ...
                                                     0.0, ...
                                                     0.0 );

            % draw to diel
            diel( :, z_coords_eme >= cur_z & z_coords_eme < cur_z + z_in ) = gratingcell_in.N;
            % update z
            cur_z               = cur_z + z_in;
            [~, cur_z_indx ]    = min( abs( z_coords_eme - cur_z ) );   % convert to array index
            
            % draw each cell
            for ii = 1:n_cells
               
                % TRICK - i can use the twoLeveLgratingcell to build the
                % dielectric for the emeSim
                gratingcell = obj.h_makeGratingCell(  obj_as_struct, ...
                                                     periods_synth(ii), ...
                                                     fills_synth(ii), ...
                                                     ratios_synth(ii), ...
                                                     offsets_synth(ii) );                                                

                % draw to diel
                try
                    diel( :, cur_z_indx:( cur_z_indx + size(gratingcell.N,2) - 1) ) = gratingcell.N;
                catch ME
                    fprintf('ERROR dielectric sizes don''t match ya better debug this bud\n');
                    error(ME);
                end
                % update z
                cur_z_indx = cur_z_indx + size(gratingcell.N,2);
                
            end
            
            % replace the dielectric in the eme object
            eme_obj.diel = diel;
            
%             % DEBUG plot the diel
%             eme_obj.plotDiel();
            
            % run EME sim
            % Converts the dielectric distribution into layers for eigen mode expansion
            eme_obj = eme_obj.convertDiel();   
            % Runs simulation
            eme_obj = eme_obj.runSimulation('plotSource','yes');      
            % compute fiber overlap
            eme_obj = eme_obj.fiberOverlap( 'zOffset', 0:.1:12,...
                                            'angleVec', -45:1:45,...
                                            'MFD', MFD * obj.units.scale * um,...
                                            'overlapDir', obj.coupling_direction);
                                        
            % DEBUG show results
            gratingUI(eme_obj);
            
            % save the final design
            final_design.directivities      = max_directivities_synth;
            final_design.fills              = fills_synth;
            final_design.ratios             = ratios_synth;
            final_design.offsets            = offsets_synth;
            final_design.periods            = periods_synth;
            final_design.angles             = angles_synth;
            final_design.scatter_strengths  = scatter_strengths_synth;
            final_design.max_coupling_eff   = eme_obj.fiberCoup.optCoup;
            final_design.MFD                = MFD;                              % in units of 'units'
            final_design.desired_angle      = angle;
            final_design.max_coupling_angle = eme_obj.fiberCoup.optAngle;
            final_design.max_coupling_offset = eme_obj.fiberCoup.optZOffset;
            final_design.power_reflection   = eme_obj.scatterProperties.PowerRefl(1,1);
            final_design.eme_obj            = eme_obj;
            obj.final_design                = final_design;
            
            
        end     % end synthesizeGaussianGrating()
        
        
        function obj  = sweepPeriodFill(obj)
            % DEPRECATED
            % Sweeps extrema of period and fill to get sense of possible
            % angular distribution of these grating cell dimensions
            %
            % 4 cases: min period, one layer
            %          max period, one layer
            %          min period, two layers
            %          max period, two layers
            
            tic;
            
            % grab parameters
            min_period  = min( obj.period_vec(:) );
            max_period  = max( obj.period_vec(:) );
%             fills       = obj.fill_vec;             % maybe change this
            fills = linspace(0.1, 0.9, 18);
%             fills = 0.5;
            
%             % sweep min period, two layers
%             fprintf('Sweep 1 of 4\n\n');
%             dir_min_two     = zeros( size(fills) );     % directivities
%             scatter_min_two = zeros( size(fills) );     % scatter strengths
%             angles_min_two  = zeros( size(fills) );     % angles
%             for ii = 1:length(fills)
%                
%                 % make grating cell
%                 GC = makeGratingCell( obj, min_period, fills(ii), 1.0, 0 );
%                 
%                 % run simulation
%                 GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
%                 
%                 % save parameters
%                 if strcmp(obj.coupling_direction, 'up')
%                     % coupling direction is up
%                     dir_min_two(ii)     = GC.directivity;
%                     scatter_min_two(ii) = GC.alpha_up;
%                     angles_min_two(ii)  = GC.max_angle_up;
%                 else
%                     % coupling direction is down
%                     dir_min_two(ii)     = 1/GC.directivity;
%                     scatter_min_two(ii) = GC.alpha_down;
%                     angles_min_two(ii)  = GC.max_angle_down;
%                 end
%                 
%                 toc;
%             end
%             
%             % DEBUG
%             GC.plotIndex();
%             GC.k
%             
%             % sweep max period, two layers
%             fprintf('Sweep 2 of 4\n\n');
%             dir_max_two     = zeros( size(fills) );     % directivities
%             scatter_max_two = zeros( size(fills) );     % scatter strengths
%             angles_max_two  = zeros( size(fills) );     % angles
%             for ii = 1:length(fills)
%                
%                 % make grating cell
%                 GC = makeGratingCell( obj, max_period, fills(ii), 1.0, 0 );
%                 
%                 % run simulation
%                 GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
%                 
%                 % save parameters
%                 if strcmp(obj.coupling_direction, 'up')
%                     % coupling direction is up
%                     dir_max_two(ii)     = GC.directivity;
%                     scatter_max_two(ii) = GC.alpha_up;
%                     angles_max_two(ii)  = GC.max_angle_up;
%                 else
%                     % coupling direction is down
%                     dir_max_two(ii)     = 1/GC.directivity;
%                     scatter_max_two(ii) = GC.alpha_down;
%                     angles_max_two(ii)  = GC.max_angle_down;
%                 end
%                 
%                 toc;
%             end
% %             
%             % sweep min period, one layers
%             fprintf('Sweep 3 of 4\n\n');
%             dir_min_one     = zeros( size(fills) );     % directivities
%             scatter_min_one = zeros( size(fills) );     % scatter strengths
%             angles_min_one  = zeros( size(fills) );     % angles
%             for ii = 1:length(fills)
%                
%                 % make grating cell
%                 GC = makeGratingCell( obj, min_period, fills(ii), 0, 0 );
%                 
%                 % run simulation
%                 GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
%                 
%                 % save parameters
%                 if strcmp(obj.coupling_direction, 'up')
%                     % coupling direction is up
%                     dir_min_one(ii)     = GC.directivity;
%                     scatter_min_one(ii) = GC.alpha_up;
%                     angles_min_one(ii)  = GC.max_angle_up;
%                 else
%                     % coupling direction is down
%                     dir_min_one(ii)     = 1/GC.directivity;
%                     scatter_min_one(ii) = GC.alpha_down;
%                     angles_min_one(ii)  = GC.max_angle_down;
%                 end
%                
%                 toc;
%             end
%             
%             % DEBUG
%             GC.plotIndex();
%             GC.k
%             
%             % sweep max period, one layers
%             fprintf('Sweep 4 of 4\n\n');
%             dir_max_one     = zeros( size(fills) );     % directivities
%             scatter_max_one = zeros( size(fills) );     % scatter strengths
%             angles_max_one  = zeros( size(fills) );     % angles
%             for ii = 1:length(fills)
%                
%                 % make grating cell
%                 GC = makeGratingCell( obj, max_period, fills(ii), 0, 0 );
%                 
%                 % run simulation
%                 GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
%                 
%                 % save parameters
%                 if strcmp(obj.coupling_direction, 'up')
%                     % coupling direction is up
%                     dir_max_one(ii)     = GC.directivity;
%                     scatter_max_one(ii) = GC.alpha_up;
%                     angles_max_one(ii)  = GC.max_angle_up;
%                 else
%                     % coupling direction is down
%                     dir_max_one(ii)     = 1/GC.directivity;
%                     scatter_max_one(ii) = GC.alpha_down;
%                     angles_max_one(ii)  = GC.max_angle_down;
%                 end
%                 
%                 toc;
%             end
%             
%             % plot angles for min period, two layers
%             figure;
%             plot( fills, angles_min_two, '-o' );
%             xlabel('fill'); ylabel('angle (deg)');
%             title(['Angles vs. fill, min period ', num2str(min_period), ' two levels']);
%             makeFigureNice();
%             
%             % plot angles for max period, two layers
%             figure;
%             plot( fills, angles_max_two, '-o' );
%             xlabel('fill'); ylabel('angle (deg)');
%             title(['Angles vs. fill, max period ', num2str(max_period), ' two levels']);
%             makeFigureNice();
%             
%             % plot angles for min period, one layers
%             figure;
%             plot( fills, angles_min_one, '-o' );
%             xlabel('fill'); ylabel('angle (deg)');
%             title(['Angles vs. fill, min period ', num2str(min_period), ' one level']);
%             makeFigureNice();
%             
%             % plot angles for max period, one layers
%             figure;
%             plot( fills, angles_max_one, '-o' );
%             xlabel('fill'); ylabel('angle (deg)');
%             title(['Angles vs. fill, max period ', num2str(max_period), ' one levels']);
%             makeFigureNice();
            
            % Sweep periods instead
            fprintf('Sweep period\n\n');
            fill                = 0.2;
            periods             = 500:20:1300;
            dir_v_period        = zeros( size(periods) );     % directivities
            scatter_v_period    = zeros( size(periods) );     % scatter strengths
            angles_v_period     = zeros( size(periods) );     % angles
            for ii = 1:length(periods)
               
                fprintf('loop %i of %i\n', ii, length(periods));
                
                % make grating cell
                GC = makeGratingCell( obj, periods(ii), fill, 1.0, 0 );
                
                % run simulation
                GC = GC.runSimulation( obj.modesolver_opts.num_modes, obj.modesolver_opts.BC, obj.modesolver_opts.pml_options );
                
                % save parameters
                if strcmp(obj.coupling_direction, 'up')
                    % coupling direction is up
                    dir_v_period(ii)        = GC.directivity;
                    scatter_v_period(ii)    = GC.alpha_up;
                    angles_v_period(ii)     = GC.max_angle_up;
                else
                    % coupling direction is down
                    dir_v_period(ii)        = 1/GC.directivity;
                    scatter_v_period(ii)    = GC.alpha_down;
                    angles_v_period(ii)     = GC.max_angle_down;
                end
                
                toc;
            end
            
            % plot angles vs. period
            figure;
            plot( periods, angles_v_period, '-o' );
            xlabel('period'); ylabel('angle (deg)');
            title(['Angles vs. period, two levels, fill ' num2str(fill)]);
            makeFigureNice();
            % plot scatter vs. period
            figure;
            plot( periods, scatter_v_period, '-o' );
            xlabel('period'); ylabel('scatter strength');
            title(['Scatter strength vs. period, two levels, fill ' num2str(fill)]);
            makeFigureNice();
            
        end     % end sweepPeriodFill()
        
        
        function [obj, GC] = testMakeGratingCell( obj, period, fill, ratio, offset_ratio )
            % TEST/DEBUGGING function for testing the drawing of a grating cell 
            % also runs the simulation lol
            %
            % this function is somewhat deprecated as of 12/7/17
            
            % make grating cell
            GC = makeGratingCell( obj, period, fill, ratio, offset_ratio );
            
            % run bloch complex k modesolver and return values
            num_modes   = 20;
            BC          = 0;                    % 0 for PEC, 1 for PMC
            pml_options = [ 1, 200, 500, 2 ];   %  [ yes/no, length in nm, strength, pml poly order ]
            
            % run simulation
            GC = GC.runSimulation( num_modes, BC, pml_options );
            
%             % plot stuff
%             GC.plotIndex();
%             GC.plotEz();
        end
        

    end     % End methods section
    
end     % end class definition



% -------------------------------------------------------------------------
% Begin auxiliary non-class methods
% -------------------------------------------------------------------------

function GC = makeGratingCell( synth_obj, period, fill_top, fill_bot, offset_ratio )
            % makes and returns a c_twoLevelGratingCell object
            % 
            % inputs:
            %   synth_obj
            %       type: c_synthGrating object AS STRUCT
            %       desc: c_synthGrating object AS STRUCT
            %   period
            %       type: double, scalar
            %       desc: period of the grating cell
            %   fill - OLD
            %       type: double, scalar
            %       desc: ratio of bottom layer to period
            %   ratio - OLD
            %       type: double, scalar
            %       desc: ratio of top layer to bottom layer
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
            
            
            % set domain 
            domain_size     = synth_obj.domain_size;
            domain_size(2)  = period;
            
            % make grating cell
            GC = c_twoLevelGratingCell( 'discretization', synth_obj.discretization, ...
                                        'units', synth_obj.units.name, ...
                                        'lambda', synth_obj.lambda, ...
                                        'domain_size', domain_size, ...
                                        'background_index', synth_obj.background_index );

            % draw cell
            % draw two levels using two level builder function
            % the inputs are organized [ top level, bottom level ]
            wg_thick        = synth_obj.waveguide_thicks;
            wg_min_y        = [ domain_size(1)/2, domain_size(1)/2-wg_thick(1) ];
%             wgs_duty_cycles = [ fill*ratio, fill ];
            wgs_duty_cycles = [ fill_top, fill_bot ];
            wgs_offsets     = [ 0, offset_ratio*period ];
            GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, synth_obj.waveguide_index, ...
                                                    wgs_duty_cycles, wgs_offsets );
            
end




















































