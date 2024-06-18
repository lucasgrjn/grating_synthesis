classdef c_synthHighDirectivityGrating
% Synthesizes a 2-level grating in an arbitrary process
% save state: 6/11/18
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
%   The user should define their own custom grating unit cell
%   drawing function.
%   HOWEVER, this function MUST have the following inputs and outputs, IN
%   ORDER:
%       function GC = your_makeGratingCell_function( synth_obj, period, fill_top, fill_bot, offset_ratio )
%           % makes and returns a c_twoLevelGratingCell object
% 
%           inputs:
%               synth_obj
%               type: c_synthGrating object AS STRUCT
%               desc: c_synthGrating object AS STRUCT
%           period
%               type: double, scalar
%               desc: period of the grating cell
%           fill_top
%               type: double, scalar
%               desc: ratio of top layer to period
%           fill_bot
%               type: double, scalar
%               desc: ratio of bottom layer to period
%           offset_ratio
%               type: double, scalar
%               desc: ratio of bottom layer offset to period
%       outputs:
%           GC
%               type: c_twoLevelGratingCell object
%               desc: two level grating cell object
%
%   As an example:
%       see the default makeGratingCell() function that's included at the
%       end of this class definition. It makes a simple waveguide + air
%       cladding grating, no fancyness
%
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
%   'optimal_angle'
%       type: double, scalar
%       desc: desired output angle, in deg
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

        discretization;     % dx and dy
        units;              % units, verbose, 'm' or 'mm', or 'um', or 'nm'
                            % has fields 'name' and 'scale'
        lambda;             % center wavelength
        background_index;   % background index
        domain_size;        % domain size, [ y size (height), x size (length) ]
        inputs;             % saves input settings for user reference

        start_time;         % time when object was created, 'YEAR-month-day hour-min-sec'
        
        coupling_direction; % direction of coupling, either 'up', or 'down'
                            % defaults to 'down'
                            
        data_directory;     % path of data directory
        data_filename;      % name of data file
        data_notes;         % verbose notes of what current sweep is doing
        data_mode;          % flag telling whether grating is run from scratch or run from previous data
                            % either 'new' or 'load'
        
        u;                  % gaussian profile, not sure if this will stay a property tho
        
        num_par_workers;    % number of parallel workers to use THIS IS DEPRECATED
        
%         modesolver_opts;    % STRUCT that stores the modesolver options
                            % CURRENTLY hardcoded.
                            
        final_design;       % STRUCT that stores the final design parameters
                            % as well as performance values
                            
        h_makeGratingCell;  % handle to the grating cell making function
        
        % parameters to optimize for
        optimal_angle;      % angle to optimize for, deviation from the normal, in deg.
        input_wg_type;      % type of input waveguide, currently supports 'bottom', 'full'
        
        % temporarily? saving the resulting variables from synthesis
        directivities_vs_fills 
        angles_vs_fills       
        periods_vs_fills        
        offsets_vs_fills        
        scatter_str_vs_fills
        k_vs_fills     
        GC_vs_fills     % this variable should be temporary, because it takes up a ton of memory
        fill_tops
        fill_bots;
        offsets;
        fill_top_bot_ratio;
        dir_b4_period_vs_fills;
        
        % final synth results
        dir_synth
        bot_fill_synth
        top_bot_fill_ratio_synth
        period_synth
        offset_synth
        angles_synth
        scatter_str_synth
        k_synth
        GC_synth    
        des_scatter_synth
        final_index                 % final index distribution
        
        
        % struct that holds debug field
        % currently has fields: final_index
%         debug;

                            
    end
    
    methods
        
        function obj = c_synthGrating(varargin)
            % Constructor
            % See top comments for input documentation
            
            % Dependency imports
            fname           = mfilename;                                            % name of class
            fpath           = mfilename('fullpath');                                % full path, including fname
            projectpath     = erase( fpath, [ 'main' filesep fname] );              % now only holds path to project's code
            % path to emeSim
            addpath([ projectpath 'eme' ]);
            % path to parfor progress monitor
            addpath([ projectpath 'auxiliary_functions' filesep 'ParforProgMon' ]);
            

            % inputs and defaults
            inputs = {  'discretization',   'none', ...
                        'units',            'nm',   ...
                        'lambda',           'none', ...
                        'background_index', 1.0,    ...
                        'domain_size',      'none', ...
                        'optimal_angle',    'none', ...
                        'coupling_direction', 'down', ...
                        'data_directory',   '', ...
                        'data_filename',    '', ...
                        'data_notes',       '', ...
                        'data_mode',        'new', ...
                        'num_par_workers',  'none', ...
                        'h_makeGratingCell', @makeGratingCell ...
                     }; 
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
                obj.start_time = datestr( datetime('now'), 'yyyy_mm_dd_HH_MM_SS_' );

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

                % set other properties
                obj.discretization      = p.discretization;
                obj.lambda              = p.lambda;
                obj.background_index    = p.background_index;
                obj.domain_size         = p.domain_size;
                obj.optimal_angle       = p.optimal_angle;

                if strcmp( p.coupling_direction, 'up') || strcmp( p.coupling_direction, 'down') 
                    % set coupling direction
                    obj.coupling_direction = p.coupling_direction;
                else
                    error('Error: input ''coupling_direction'' is not valid. Valid entries are ''up'' or ''down''. You entered ''%s''', p.coupling_direction);
                end

                % set file saving/loading properties
                obj.data_directory  = p.data_directory;
                obj.data_filename   = p.data_filename;
                obj.data_notes      = p.data_notes;
                obj.data_mode       = p.data_mode;
                
                % number of parallel workers
                obj.num_par_workers = p.num_par_workers;
                
%                 % default modesolver options (currently hardcoded)
%                 num_modes   = 20;
%                 BC          = 0;                    % 0 for PEC, 1 for PMC
%                 pml_options = [ 1, 200, 500, 2 ];   % [ yes/no, length in nm, strength, pml poly order ]
%                 obj.modesolver_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options );
                
                % set handle to grating cell making function
                obj.h_makeGratingCell = p.h_makeGratingCell;
                
                
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
            % This function will be DEPRECATED
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
            % This function will either be deprecated or will need to be
            % updated to work with the newer synthesis pipeline
            % 
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
        
        
        
        function obj = synthesizeUniformGrating(obj, MFD, fill_factor_top, fill_factor_bot, input_wg_type, DEBUG)
            % Synthesizes a uniform grating at the desired angle
            %
            %   User supplies the waveguide fill/period ratio. The fill
            %   ratios determine the scattering strength/reflection of the
            %   final grating
            %
            % Inputs:
            %   MFD
            %       type: scalar, double
            %       desc: mode field diameter of coupling fiber, in units 'units'
            %   fill_factor_top
            %       type: scalar, double
            %       desc: top fill factor (ratio of silicon/period)
            %   fill_factor_bot
            %       type: scalar, double
            %       desc: bottom fill factor (ratio of silicon/period)
            %   input_wg_type
            %       type: string
            %       desc: 'normal' for body + poly input waveguide
            %             'invert' for body input waveguide only
            %             This is kind of an odd way to specify the input
            %             waveguide type... I'll figure out a cleaner/more
            %             intuitive convention later
            %   DEBUG
            %       type: boolean
            %       desc: OPTIONAL flag - set to true to enable debug mode
            
            
            % default debug mode to false
            if nargin < 7
                DEBUG = false;
            end
           
            % get waveguide k
            fprintf('Simulating waveguide...\n');         

            % make grating cell
            waveguide = obj.h_makeGratingCell( obj.convertObjToStruct(), obj.discretization, 1.0, 1.0, 0.0 );
            
            % run simulation
            % sim settings
            lambda_nm   = obj.lambda * obj.units.scale * 1e9;                               % units nm
            guess_n     = 0.7 * max( waveguide.N(:) );                                      % guess index
            guessk      = guess_n * 2*pi/lambda_nm;                                         % units 1/nm
            num_modes   = 5;
            BC          = 0;                                                                % 0 = PEC
            pml_options = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
            % run sim
            waveguide   = waveguide.runSimulation( num_modes, BC, pml_options, guessk );
            
            % update guessk (units 1/nm)
            guessk = waveguide.k;
            
            % grab waveguide k
            waveguide_k = waveguide.k * obj.units.scale * 1e9;                              % in units 1/'units'                          
            
%             % DEBUG plot stuff
%             waveguide.plotEz_w_edges();
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0      = obj.background_index * ( 2*pi/obj.lambda );
            kx      = k0 * sin( (pi/180) * obj.optimal_angle );
            period  = 2*pi/(waveguide_k- kx);                                               % units of 'units'
            
            % snap units to discretization
            guess_period    = obj.discretization * round(period/obj.discretization);
            guess_period_nm = guess_period * obj.units.scale * 1e9;
            disc_nm         = obj.discretization * obj.units.scale * 1e9;                   % discretization in nm
            fprintf('...done\n\n');
            
            % sweep FF
            fprintf('Sweeping fill factors...\n');
            
            % set fill factors
            step_size   = 0.025;                                            % amount to decrease fill factor by each time (approx)
            max_fill    = 0.95;                                             % max fill
            n_fills     = round(max_fill/step_size);                        % number of fill factor combinations to try
            fill_tops   = linspace( max_fill, fill_factor_top, n_fills );
            fill_bots   = linspace( max_fill, fill_factor_bot, n_fills );
            
            
            % initialize saving variables
            % not sure if each of these will be used, but they are at least
            % useful for debugging purposes
            directivities_vs_fills  = zeros( 1, n_fills );    
            angles_vs_fills         = zeros( 1, n_fills );     
            periods_vs_fills        = zeros( 1, n_fills ); 
            offsets_vs_fills        = zeros( 1, n_fills ); 
            scatter_str_vs_fills    = zeros( 1, n_fills ); 
            k_vs_fills              = zeros( 1, n_fills ); 
            GC_vs_fills             = cell( 1, n_fills ); 
            
            % set solver settings
            num_modes   = 1;
            BC          = 0;                                                % 0 = PEC
            pml_options = [1, 200, 20, 2]; 
            
            tic;
            
            % sweep, optimize the period and guessk for offset = 0
            for i_fill = 1:n_fills
                % for each top/bottom fill combination
                    
                % print iteration
                fprintf('Fill factor iteration %i of %i\n', i_fill, n_fills );

                % sweep periods, in nm
                % only sweep larger periods. Doubtful that the period
                % will be smaller
                periods_nm     = guess_period_nm : disc_nm : 1.1 * guess_period_nm;
                periods_nm     = disc_nm * round( periods_nm/disc_nm );

                % init saving variables
                angles          = zeros( size(periods_nm) );
                k_vs_period     = zeros( size(periods_nm) );
                GC_vs_period    = cell( size(periods_nm) );

                % sweep periods
                fprintf('Sweeping periods...\n');
                for i_period = 1:length(periods_nm)

                    fprintf('Iteration %i of %i\n', i_period, length(periods_nm) );

                    % make grating cell
                    GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                                periods_nm(i_period), ...
                                                fill_tops(i_fill), ...
                                                fill_bots(i_fill), ...
                                                0.0 );

                    % run sim
                    GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

                    % save angle
                    if strcmp( obj.coupling_direction, 'up' )
                        % coupling direction is upwards
                        angles( i_period ) = GC.max_angle_up;
                    else
                        % coupling direction is downwards
                        angles( i_period ) = GC.max_angle_down;
                    end

                    % update GC list
                    GC_vs_period{i_period} = GC;

                    % update k
                    k_vs_period(i_period)   = GC.k;
                    guessk                  = GC.k;

                    toc;

                end
                fprintf('...done.\n');


                % pick best period
                [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles ) );
                best_period_nm                  = periods_nm( indx_best_period );
                best_period_k                   = k_vs_period( indx_best_period );
                best_GC                         = GC_vs_period{ indx_best_period };


                % save data
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    directivities_vs_fills( i_fill )   = best_GC.directivity;
                    angles_vs_fills( i_fill )          = best_GC.max_angle_up;
                    scatter_str_vs_fills( i_fill )     = best_GC.alpha_up;
                else
                    % coupling direction is downwards
                    directivities_vs_fills( i_fill )   = 1./best_GC.directivity;
                    angles_vs_fills( i_fill )          = best_GC.max_angle_down;
                    scatter_str_vs_fills( i_fill )     = best_GC.alpha_down;
                end
                periods_vs_fills( i_fill )  = best_period_nm * 1e-9 / obj.units.scale;  % in units 'units'
                k_vs_fills( i_fill )        = best_GC.k * 1e9 * obj.units.scale;        % in units 1/'units'
                GC_vs_fills{ i_fill }       = best_GC;


                % update the period and the guessk
                guessk              = best_GC.k;
                guess_period_nm     = best_period_nm;

                
            end     % end for ii = 1:n_fills
            fprintf('..done sweeping fill factors\n\n');
            
            % save variables to object
            % optional for this method
            obj.directivities_vs_fills  = directivities_vs_fills;
            obj.angles_vs_fills         = angles_vs_fills;
            obj.scatter_str_vs_fills    = scatter_str_vs_fills;
            obj.periods_vs_fills        = periods_vs_fills;
            obj.offsets_vs_fills        = offsets_vs_fills;
            obj.k_vs_fills              = k_vs_fills;
            obj.GC_vs_fills             = GC_vs_fills;
            
            % optimize offset
            
            % range of offsets to try
            offsets = 0:0.02:0.98;
            
            % init vars
            directivities = zeros(1, length(offsets));                      % dimensions vs. offsets
            GC_vs_offsets = cell(1, length(offsets));
            
            % sweep offsets
            fprintf('Sweeping offsets...\n');
            for i_offset = 1:length(offsets)
               
                fprintf('Iteration %i of %i\n', i_offset, length(offsets) );

                    % make grating cell
                    GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                                guess_period_nm, ...
                                                fill_factor_top, ...
                                                fill_factor_bot, ...
                                                offsets(i_offset) );

                    % run sim
                    GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
                    
                    % save directivity
                    if strcmp( obj.coupling_direction, 'up' )
                        % coupling direction is upwards
                        directivities( i_offset ) = GC.directivity;
                    else
                        % coupling direction is downwards
                        directivities( i_offset ) = 1./( GC.directivity );
                    end

                    % update variables
                    guessk                      = GC.k;
                    GC_vs_offsets{ i_offset }   = GC;

                    toc;
                
            end     % end sweeping offsets
            fprintf('...done\n\n');
            
            % pick best offset
            [ ~, indx_best_offset ]     = max( directivities );
            best_offset                 = offsets( indx_best_offset );
            best_GC                     = GC_vs_offsets{ indx_best_offset };
            
            
            % run one more period optimization
            
            % grab angle
            if strcmp( obj.coupling_direction, 'up' )
                % coupling direction is upwards
                angle_after_offset = best_GC.max_angle_up;
            else
                % coupling direction is downwards
                angle_after_offset = best_GC.max_angle_down;
            end
            
            % depending on angle, sweep larger or smaller periods
            if angle_after_offset > obj.optimal_angle
                % phase matching too short, sweep shorter periods
                periods_nm_sweep2 = guess_period_nm : -disc_nm : 0.8 * guess_period_nm;
            else
                % phase matching too long, sweep longer periods
                periods_nm_sweep2 = guess_period_nm : disc_nm : 1.2 * guess_period_nm;
            end
            % snap to grid
            periods_nm_sweep2 = disc_nm * round( periods_nm_sweep2/disc_nm );
            
            % reset angles
            angles_sweep2       = zeros( 1, length(periods_nm_sweep2) );
            GC_vs_period_sweep2 = cell( 1, length(periods_nm_sweep2) );
            
            
            % sweep periods
            fprintf('Sweeping periods again...\n');
            for i_period = 1:length(periods_nm_sweep2)

                fprintf('Iteration %i of %i\n', i_period, length(periods_nm_sweep2) );

                % make grating cell
                GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                            periods_nm_sweep2(i_period), ...
                                            fill_tops(i_fill), ...
                                            fill_bots(i_fill), ...
                                            best_offset );

                % run sim
                GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

                % save angle
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    angles_sweep2( i_period ) = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    angles_sweep2( i_period ) = GC.max_angle_down;
                end

                % update GC list
                GC_vs_period_sweep2{i_period} = GC;

                % update k
                guessk = GC.k;

                toc;

            end
            fprintf('...done.\n\n');
            
            % pick best period
            [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles_sweep2 ) );
            best_period_nm                  = periods_nm_sweep2( indx_best_period );
            best_GC_sweep2                  = GC_vs_period_sweep2{ indx_best_period };
            
            
            % finally, fine tune with a local optimization
            
            
            % inputs to merit function
            weights         = [10, 1];
            fill_factors    = [ fill_factor_top, fill_factor_bot ];
            period          = best_period_nm * 1e-9 / obj.units.scale;
            guessk_nm       = guessk;
            guessk          = best_GC_sweep2.k * 1e9 * obj.units.scale;

%             % DEBUG run FOM
%             [ FOM ] = obj.merit_period_offset( [1, best_offset], weights, fill_factors, period, guessk );
            
            % starting point
            x0 = [ 1, best_offset ];

            % options
            if DEBUG
                plot_functions = {@optimplotfval, @optimplotx};
            else
                plot_functions = {};
            end
            opts = optimset( 'Display', 'iter', ...
                             'FunValCheck', 'off', ...
                             'MaxFunEvals', 400, ...
                             'MaxIter', 400, ...
                             'PlotFcns', plot_functions );
            
                         
            % run fminsearch, simplex search
            fprintf('Running local optimizer...\n');
            [x, fval, exitflag, output] = fminsearch( @(x) obj.merit_period_offset( x, weights, fill_factors, period, guessk ), x0, opts );
            toc;
            fprintf('...done\n\n');

%             % run fminsearch, simplex search
%             fprintf('Running local optimizer...\n');
%             tic;
%             [x1, fval1, exitflag, output] = fminsearch( @(x) obj.merit_period_offset( x, weights, fill_factors, period, guessk ), x0, opts );
%             toc;
%             fprintf('...done\n\n');
% 
%             % run fminunc, gradient search
%             fprintf('Running local optimizer...\n');
%             tic;
%             [x2, fval2, exitflag, output] = fminunc( @(x) obj.merit_period_offset( x, weights, fill_factors, period, guessk ), x0, opts );
%             toc;
%             fprintf('...done\n\n');
%             
%             
%             % pick the better option
%             if fval1 < fval2
%                 % simplex wins
%                 x = x1;
%             else
%                 % gradient wins
%                 x = x2;
%             end
            
            
            % make final grating cell
            best_period     = x(1) * period;
            best_period     = obj.discretization * round( best_period/obj.discretization );     % snap to grid
            best_period_nm  = best_period * obj.units.scale * 1e9;
            best_offset     = x(2);
            GC_final = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                                best_period_nm, ...
                                                fill_factor_top, ...
                                                fill_factor_bot, ...
                                                best_offset );

            % run sim
            GC_final = GC_final.runSimulation( num_modes, BC, pml_options, guessk_nm );
            
            
            % NOW to verify the design
            % Run it in EME
            fprintf('Verifying design in EME...\n');
            
            % Set Up Simulation
            % note that emeSim uses 'z' as propagation direction and 'x'
            % as transverse (synthGrating uses 'x' and 'y' respectively)
            % and units are in um
            n_cells     = 20;
            um          = 1e6;
            dx          = obj.discretization * obj.units.scale * um;                % in um
            dz          = 5e-3;                                                     % in um
            pol         = 0;                                                        % 0 for TE, 1 for TM
            z_in        = (MFD * obj.units.scale * um)/2;                           % length of input section of waveguide
            xf          = obj.domain_size(1) * obj.units.scale * um;                % in um
            zf          = n_cells * best_period * obj.units.scale * um + z_in;      % in um
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
            if strcmp( input_wg_type, 'invert' )
                % invert design, input is body wg only
                
                gratingcell_in  = obj.h_makeGratingCell( obj_as_struct, ...
                                                     z_in/(um*obj.units.scale), ...
                                                     0.0, ...
                                                     1.0, ...
                                                     0.0 );
                
            elseif strcmp( input_wg_type, 'normal' )
                % normal design, input is body + poly wg
                
                gratingcell_in  = obj.h_makeGratingCell( obj_as_struct, ...
                                                     z_in/(um*obj.units.scale), ...
                                                     1.0, ...
                                                     1.0, ...
                                                     0.0 );
                
            else
                % throw error, input_wg_type was invalid
                error('input_wg_type must either be "invert" or "normal"');
            end


            % draw to diel
            diel( :, z_coords_eme >= cur_z - dz/10 & z_coords_eme < cur_z + z_in - dz/10 ) = gratingcell_in.N;
            % update z
            cur_z               = cur_z + z_in;
            [~, cur_z_indx ]    = min( abs( z_coords_eme - cur_z ) );   % convert to array index
            
            % draw each cell
            gratingcell = obj.h_makeGratingCell(  obj_as_struct, ...
                                                     best_period, ...
                                                     fill_factor_top, ...
                                                     fill_factor_bot, ...
                                                     best_offset );
            gratingcell_index_rep = repmat( gratingcell.N, 1, n_cells );
            
            % replace the dielectric in the eme object
            diel( :, cur_z_indx:end )    = gratingcell_index_rep;
            eme_obj.diel = diel;
            
            % DEBUG plot the diel
            if DEBUG
                eme_obj.plotDiel();
            end
            
            % run EME sim
            % Converts the dielectric distribution into layers for eigen mode expansion
            eme_obj = eme_obj.convertDiel();   
            % Runs simulation
            eme_obj = eme_obj.runSimulation('plotSource','no');      
            % compute fiber overlap
            eme_obj = eme_obj.fiberOverlap( 'zOffset', 0:.1:20,...
                                            'angleVec', -45:1:45,...
                                            'MFD', MFD * obj.units.scale * um,...
                                            'overlapDir', obj.coupling_direction);
                                        
%             % DEBUG show results
%             gratingUI(eme_obj);
            
            % store results
            obj.final_design.GC_final               = GC_final;
            obj.final_design.eme_obj                = eme_obj;
            obj.final_design.max_coupling_eff       = eme_obj.fiberCoup.optCoup;
            obj.final_design.max_coupling_offset    = eme_obj.fiberCoup.optZOffset / ( um * obj.units.scale );    % in units 'units'
            obj.final_design.max_coupling_angle     = eme_obj.fiberCoup.optAngle;
            obj.final_design.reflection_coeff       = eme_obj.scatterProperties.PowerRefl(1,1);
            obj.final_design.period                 = best_period;
            obj.final_design.offset                 = best_offset;
            obj.final_design.fill_factor_top        = fill_factor_top;
            obj.final_design.fill_factor_bot        = fill_factor_bot;
            obj.final_design.desired_angle          = obj.optimal_angle;
            obj.final_design.desired_MFD            = MFD;
            

        end         % end synthesizeUniformGrating()
        
        
        function [ FOM ] = merit_period_offset( obj, inputs, weights, fill_factors, period, guessk )
        % Merit function used to optimize grating cell's period and
        % offset
        % 
        % Inputs:
        %   inputs
        %       type: 1x2 array
        %       desc: inputs to merit function, aka the data point to get
        %             FOM of
        %             currently [ scaling of period, absolute value of
        %             offset ]
        %             such that the simulated GC has period of
        %             period*scaling of period, and offset = offset
        %   weights
        %       type: double, vector
        %       desc: 1x2 array to weigh the two objectives
        %   fill_factors
        %       type: double, vector
        %       desc: [ fill factor top, fill factor bottom ]
        %   period
        %       type: double, scalar
        %       desc: period to scale, in units 'units'
        %   guessk
        %       type: double, scalar
        %       desc: guessk, in units rad/'units'

            % parse inputs
            period      = period * inputs(1);                               % units 'units'
%             period_nm   = period * obj.units.scale * 1e9;                   % units nm
            offset      = inputs(2);
            fill_top    = fill_factors(1);
            fill_bot    = fill_factors(2);
%             guessk_nm   = guessk / ( obj.units.scale * 1e9 );               % units rad/nm 

            % make grating coupler object
            GC = obj.h_makeGratingCell( obj.convertObjToStruct(), period, fill_top, fill_bot, offset );

            % simulation settings
            num_modes   = 1;
            BC          = 0;     % 0 for PEC, 1 for PMC
            % PML_options(1): PML in y direction (yes=1 or no=0)
            % PML_options(2): length of PML layer in nm
            % PML_options(3): strength of PML in the complex plane
            % PML_options(4): PML polynomial order (1, 2, 3...)
            pml_options = [ 1, 200, 20, 2 ];

            % run simulation
            GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

            % grab angle, directivity depending on up/down coupling
            if strcmp( obj.coupling_direction, 'up' )
                % coupling direction is upwards
                angle_sim   = GC.max_angle_up;
                directivity = GC.directivity;
            else
                % coupling direction is downwards
                angle_sim   = GC.max_angle_down;
                directivity = 1/GC.directivity;
            end
            
            
            % minimize the FOM
            % FOM = error in angle + directivity
            FOM =   weights(1) * abs(( obj.optimal_angle - angle_sim )/obj.optimal_angle) - ...
                    weights(2) * log10(directivity);

        end         % end merit_period_offset()
        
        
        
        function [ FOM ] = merit_offset_directivity( obj, input, fill_factors, period, guessk )
        % Merit function used to optimize grating cell's directivity vs.
        % offset
        % 
        % Inputs:
        %   input
        %       type: scalar, double
        %       desc: one input - the offset
        %   fill_factors
        %       type: double, vector
        %       desc: [ fill factor top, fill factor bottom ]
        %   period
        %       type: double, scalar
        %       desc: period to scale, in units 'units'
        %   guessk
        %       type: double, scalar
        %       desc: guessk, in units rad/'units'

            % parse inputs
            offset      = input;
            fill_top    = fill_factors(1);
            fill_bot    = fill_factors(2);

            % make grating coupler object
            GC = obj.h_makeGratingCell( obj.convertObjToStruct(), period, fill_top, fill_bot, offset );

            % simulation settings
            num_modes   = 1;
            BC          = 0;     % 0 for PEC, 1 for PMC
            % PML_options(1): PML in y direction (yes=1 or no=0)
            % PML_options(2): length of PML layer in nm
            % PML_options(3): strength of PML in the complex plane
            % PML_options(4): PML polynomial order (1, 2, 3...)
            pml_options = [ 1, 200, 20, 2 ];

            % run simulation
            GC = GC.runSimulation( num_modes, BC, pml_options, guessk );

            % grab angle, directivity depending on up/down coupling
            if strcmp( obj.coupling_direction, 'up' )
                % coupling direction is upwards
                directivity = GC.directivity;
            else
                % coupling direction is downwards
                directivity = 1/GC.directivity;
            end
            
            
            % minimize the FOM
            % maximize directivity
            FOM =  - log10(directivity);

        end         % end merit_offset_directivity()


        
        
        function obj = synthesizeGaussianGrating_old(obj, angle, MFD)
            % DEPRECATED VERSION
            %
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
            
            
        end     % end synthesizeGaussianGrating_old()
        
        
        function obj = synthesizeGaussianGrating(obj, MFD, DEBUG)
            % Synthesizes a grating that is mode-matched to fiber gaussian
            % mode
            %
            % Coupling angle is determined by obj.optimal_angle
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter
            %             
            %             outline
            %             1. sim wg, no perturbations
            %             2. for smallest perturbation (largest ff combo)
            %                 a. calculate analytical period
            %                 b. sweep offset, using anlaytical period
            %                 c. sweep period, using best offset
            %                 d. save
            %   DEBUG
            %       type: boolean
            %       desc: OPTIONAL, set to true to enable debug mode
            %             which plots stuff and saves stuff, etc.
            
            if nargin < 3
                DEBUG = false;
            end
            
            % generate x coordinates for the gaussian mode
            % must be large enough to fit all cells + mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiberModeGaussian(    w0, zvec, xvec,...
                                                    obj.optimal_angle, d0, obj.background_index );
            
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
            
            % -------------------------------------------------------------
            % Simulation time
            
            
           
            % get waveguide k
            fprintf('Simulating waveguide...\n');
            
%             % grab modesolver options
%             % actually i'm not using this...
%             num_modes   = obj.modesolver_opts.num_modes;
%             BC          = obj.modesolver_opts.BC;        
%             pml_options = obj.modesolver_opts.pml_options;
                

%             % make grating cell
%             waveguide = obj.h_makeGratingCell( obj.convertObjToStruct(), obj.discretization, 1.0, 1.0, 0.0 );
%             
%             % run simulation
%             % sim settings
% %             lambda_nm   = obj.lambda * obj.units.scale * 1e9;                             % units nm
%             guess_n     = 0.7 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
%             guessk      = guess_n * 2*pi/obj.lambda;                                        % units rad/'units'
%             num_modes   = 5;
%             BC          = 0;                                                                % 0 = PEC
%             pml_options = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
%             % run sim
%             waveguide   = waveguide.runSimulation( num_modes, BC, pml_options, guessk );
%             
%             % update guessk (units rad/'units')
%             guessk = waveguide.k;
%             
%             % grab waveguide k
%             waveguide_k = waveguide.k;  % units of rad/'units'          % * obj.units.scale * 1e9;                              % in units 1/'units'                          
%             
%             % DEBUG plot stuff
%             waveguide.plotEz_w_edges();
%             
%             % calculate analytical period which would approximately phase
%             % match to desired output angle
%             k0      = obj.background_index * ( 2*pi/obj.lambda );
%             kx      = k0 * sin( (pi/180) * obj.optimal_angle );
%             period  = 2*pi/(waveguide_k- kx);                                               % units of 'units'
%             
%             % snap period to discretization
%             guess_period    = obj.discretization * round(period/obj.discretization);
% %             guess_period_nm = guess_period * obj.units.scale * 1e9; 
%             fprintf('...done\n\n');
            

            % ----------------------------------------------------------
%             % sweep FF TOP AND BOT
%             fprintf('Sweeping fill factors for directivity and angle...\n');
%             
%             % set fill factors and offsets
% %             fill_tops       = fliplr( 0.3:0.025:0.95 );
% %             fill_bots       = fliplr( 0.3:0.025:0.95 );
%             % DEBUG
%             fill_tops       = fliplr( 0.0:0.025:1 );
%             fill_bots       = fliplr( 0.0:0.025:1 );
%             offsets         = 0:0.04:0.98;
%             offsets_orig    = offsets;
%             
%             % save fills and offsets
%             obj.fill_tops   = fill_tops;
%             obj.fill_bots   = fill_bots;
%             obj.offsets     = offsets;
%             
%             % initialize saving variables
%             directivities_vs_fills  = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             angles_vs_fills         = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             periods_vs_fills        = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             offsets_vs_fills        = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             scatter_str_vs_fills    = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             k_vs_fills              = zeros( length( fill_tops ), length( fill_bots ) );     % dimensions top fill vs. bot fill
%             GC_vs_fills             = cell( length( fill_tops ), length( fill_bots ) );      % dimensions top fill vs. bot fill
%             
%             % set solver settings
%             num_modes   = 1;
%             BC          = 0;                                                % 0 = PEC
%             pml_options = [1, 200, 20, 2]; 
%             
%             tic;
%             ii = 0;
%             
%             % sweep
%             for i_ff_top = 1:length( fill_tops )
%                 % For each top fill factor
%                 
%                 for i_ff_bot = 1:length( fill_bots )
%                     % for each bottom fill factor
%                     
%                     % print iteration
%                     ii = ii + 1;
%                     fprintf('Fill factor iteration %i of %i\n', ii, length( fill_tops ) * length( fill_bots ) );
%                     
%                     
% %                     if ii == 1
%                         % First go-around, sweep offset and period
%                     
%                         % init saving variables
%                         directivities = zeros( size(offsets) );
%                         k_vs_offset   = zeros( size(offsets) );
%                         angles        = zeros( size(offsets) );     % for debugging
% 
%                         
% %                         % add a little bit to the guess period to avoid the
% %                         % bandgap
% %                         guess_period = guess_period * 1.025;
% %                         guess_period = obj.discretization * round(guess_period/obj.discretization);
% 
%                         % Sweep offsets, pick offset with best directivity
%                         fprintf('Sweeping offsets...\n');
%                         for i_offset = 1:length( offsets )
% 
%                             fprintf('Iteration %i of %i\n', i_offset, length(offsets) );
% 
%                             % make grating cell
%                             GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                         guess_period, ...
%                                                         fill_tops(i_ff_top), ...
%                                                         fill_bots(i_ff_bot), ...
%                                                         offsets(i_offset) );
% 
%                             % run sim
%                             GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
%                             
%                             
% %                             % DEBUG plot field at chosen iteration
% %                             if DEBUG
% %                                 if i_offset == 20
% %                                     GC.plotEz_w_edges();
% %                                 end
% %                             end
% 
%                             % save directivity
%                             if strcmp( obj.coupling_direction, 'up' )
%                                 % coupling direction is upwards
%                                 directivities( i_offset )   = GC.directivity;
%                                 angles( i_offset )          = GC.max_angle_up;
%                             else
%                                 % coupling direction is downwards
%                                 directivities( i_offset )   = 1./( GC.directivity );
%                                 angles( i_offset )          = GC.max_angle_down;
%                             end
% 
%                             % update the guessk (units rad/'units')
%                             guessk                  = GC.k;
%                             k_vs_offset( i_offset ) = GC.k;
% 
%                             toc;
% 
%                         end     % end for i_offset = ...
%                         fprintf('...done.\n');
%                         
% %                         % DEBUG plot directivity vs. offset
%                         if DEBUG
%                             figure;
%                             plot( offsets, directivities, '-o' );
%                             xlabel('offsets'); ylabel('directivities');
%                             title('DEBUG directivities vs offsets');
%                             makeFigureNice();
% %                             
% %                             figure;
% %                             plot( offsets, angles, '-o' );
% %                             xlabel('offsets'); ylabel('angles');
% %                             title('DEBUG angles vs offsets for first run');
% %                             makeFigureNice();
% %                             
%                         end
% 
%                         % pick best offset
%                         [ ~, indx_best_offset ]     = max( directivities );
%                         best_offset                 = offsets( indx_best_offset );
%                         best_offset_k               = k_vs_offset( indx_best_offset );
%                         
%                         % DEBUG plot grating with best directivity
%                         if DEBUG
%                            
%                             % make grating cell
%                             GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                         guess_period, ...
%                                                         fill_tops(i_ff_top), ...
%                                                         fill_bots(i_ff_bot), ...
%                                                         best_offset );
% 
%                             % run sim
%                             GC = GC.runSimulation( num_modes, BC, pml_options, best_offset_k );
%                             
%                             % plot field
%                             GC.plotEz_w_edges();
%                             
%                         end
% 
%                         
%                         % now sweep periods
%                         % only sweep larger periods. Doubtful that the period
%                         % will be smaller
%                         periods     = guess_period : obj.discretization : 1.05 * guess_period;
%                         periods     = obj.discretization * round(periods/obj.discretization);
% %                         periods_nm  = periods * obj.units.scale * 1e9;                            % convert to nm
% 
%                         % init saving variables
%                         angles          = zeros( size(periods) );
%                         k_vs_period     = zeros( size(periods) );
%                         GC_vs_period    = cell( size(periods) );
% 
%                         % sweep periods
%                         guessk = best_offset_k;
%                         fprintf('Sweeping periods...\n');
%                         for i_period = 1:length(periods)
% 
%                             fprintf('Iteration %i of %i\n', i_period, length(periods) );
% 
%                             % make grating cell
%                             GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                         periods(i_period), ...
%                                                         fill_tops(i_ff_top), ...
%                                                         fill_bots(i_ff_bot), ...
%                                                         best_offset );
% 
%                             % run sim
%                             GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
%                             % save angle
%                             if strcmp( obj.coupling_direction, 'up' )
%                                 % coupling direction is upwards
%                                 angles( i_period ) = GC.max_angle_up;
%                             else
%                                 % coupling direction is downwards
%                                 angles( i_period ) = GC.max_angle_down;
%                             end
% 
%                             % update GC list
%                             GC_vs_period{i_period} = GC;
% 
%                             % update k (units of rad/'units')
%                             k_vs_period(i_period)   = GC.k;
%                             guessk                  = GC.k;
% 
%                             toc;
% 
%                         end
%                         fprintf('...done.\n');
% 
%                         % pick best period
%                         [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles ) );
%                         best_period                     = periods( indx_best_period );
%                         best_period_k                   = k_vs_period( indx_best_period );
%                         best_GC                         = GC_vs_period{ indx_best_period };
%                         
% %                     else
% %                         % Run local optimizer
% % 
% %                         % grab fill factors
% %                         fill_factors    = [ fill_tops(i_ff_top), fill_bots(i_ff_bot) ];
% %                         
% %                         % First optimize the offset
% %                         % starting point
% %                         x0 = guess_offset;
% %                         
% %                         % options
% %                         opts = optimset( 'Display', 'iter', ...
% %                                          'FunValCheck', 'off', ...
% %                                          'MaxFunEvals', 400, ...
% %                                          'MaxIter', 400 );
% % 
% %                         % run fminsearch, simplex search
% %                         % returns x = [ offset ]
% %                         fprintf('Optimizing offset for max directivity...\n');
% %                         [guess_offset, fval, exitflag, output] = fminsearch( @(x) obj.merit_offset_directivity( x, fill_factors, guess_period, guessk ), x0, opts );
% %                         toc;
% %                         fprintf('...done\n\n');
% %                         
% %                         % update guessk
% %                         % make grating coupler object
% %                         GC = obj.h_makeGratingCell( obj.convertObjToStruct(), guess_period, fill_tops(i_ff_top), fill_bots(i_ff_bot), guess_offset );
% %                         % simulation settings
% %                         num_modes   = 1;
% %                         BC          = 0;     % 0 for PEC, 1 for PMC
% %                         pml_options = [ 1, 200, 20, 2 ];
% %                         % run simulation
% %                         GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% %                         % update guessk
% %                         guessk = GC.k;
% %                         
% %                         
% %                         % Now optimize both the offset and the period
% %                         % inputs to merit function
% %                         weights         = [10, 1];                                                  % angle, offset                        
% % 
% %                         % starting point
% %                         x0 = [ 1, guess_offset ];
% % 
% %                         % run fminsearch, simplex search
% %                         % returns x = [ period ratio, offset ]
% %                         fprintf('Running local optimizer...\n');
% %                         [x, fval, exitflag, output] = fminsearch( @(x) obj.merit_period_offset( x, weights, fill_factors, guess_period, guessk ), x0, opts );
% %                         toc;
% %                         fprintf('...done\n\n');
% %                         
% %                         
% %                         % parse results
% %                         best_period = x(1) * guess_period;
% %                         best_period = obj.discretization * round( best_period/obj.discretization );     % snap to grid
% %                         best_offset = x(2);
% %                         
% %                         
% %                         % finally resimulate the GC with the chosen offset and
% %                         % period and save this data
% %                         % make grating cell
% %                         best_GC = obj.h_makeGratingCell(    obj.convertObjToStruct(), ...
% %                                                             best_period, ...
% %                                                             fill_tops(i_ff_top), ...
% %                                                             fill_bots(i_ff_bot), ...
% %                                                             best_offset );
% % 
% %                         % run sim
% %                         best_GC = best_GC.runSimulation( 1, BC, pml_options, guessk );
% %                         
% %                     
% %                     end     % end if ii == 0
%                     
%                     
%                     % save data
%                     if strcmp( obj.coupling_direction, 'up' )
%                         % coupling direction is upwards
%                         directivities_vs_fills( i_ff_top, i_ff_bot )   = best_GC.directivity;
%                         angles_vs_fills( i_ff_top, i_ff_bot )          = best_GC.max_angle_up;
%                         scatter_str_vs_fills( i_ff_top, i_ff_bot )     = best_GC.alpha_up;
%                     else
%                         % coupling direction is downwards
%                         directivities_vs_fills( i_ff_top, i_ff_bot )   = 1./best_GC.directivity;
%                         angles_vs_fills( i_ff_top, i_ff_bot )          = best_GC.max_angle_down;
%                         scatter_str_vs_fills( i_ff_top, i_ff_bot )     = best_GC.alpha_down;
%                     end
%                     periods_vs_fills( i_ff_top, i_ff_bot )  = best_period;
%                     offsets_vs_fills( i_ff_top, i_ff_bot )  = best_offset;
%                     k_vs_fills( i_ff_top, i_ff_bot )        = best_GC.k;
%                     GC_vs_fills{ i_ff_top, i_ff_bot }       = best_GC;
%                     
%                     
%                     % update the guess parameters, period, k, offset
%                     if i_ff_bot == 1
%                         % first iteration, save these guess parameters for
%                         % the next top level loop
%                         next_top_loop_period    = best_period;
%                         next_top_loop_k         = best_GC.k;
%                         next_top_loop_offset    = best_offset;
%                     end
%                     guessk              = best_GC.k;
%                     guess_period        = best_period;
%                     guess_offset        = best_offset;
%                     
%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                     
%                     
%                 end     % end for i_ff_bot = ...
%                 
%                 % update the guess parameters, period, k, offset
%                 guess_period    = next_top_loop_period;
%                 guessk          = next_top_loop_k;
%                 guess_offset    = next_top_loop_offset;
%                 
%                 % update the offsets
%                 % grab previous offset index
%                 [~, indx_prev_offset] = min( abs( offsets_orig - offsets_vs_fills( i_ff_top, 1 ) ) );
%                 % shift offsets to start at previous offset
%                 offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                 
%             end     % end for i_ff_top = ...
%             
%             % save variables to object
%             obj.directivities_vs_fills  = directivities_vs_fills;
%             obj.angles_vs_fills         = angles_vs_fills;
%             obj.scatter_str_vs_fills    = scatter_str_vs_fills;
%             obj.periods_vs_fills        = periods_vs_fills;
%             obj.offsets_vs_fills        = offsets_vs_fills;
%             obj.k_vs_fills              = k_vs_fills;
%             obj.GC_vs_fills             = GC_vs_fills;
%             
%             % DEBUG plot stuff
%             if DEBUG == true
%                 % directivity vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, 10*log10(directivities_vs_fills) );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Directivity (dB) vs. fill factors');
%                 savefig('dir_v_ff.fig');
%                 saveas(gcf, 'dir_v_ff.png');
% 
%                 % angles vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, angles_vs_fills );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Angles (deg) vs. fill factors');
%                 savefig('angle_v_ff.fig');
%                 saveas(gcf, 'angle_v_ff.png');
% 
%                 % scattering strength alpha vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, real(scatter_str_vs_fills) );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Scattering strength (real) vs. fill factors');
%                 savefig('scatter_str_v_ff.fig');
%                 saveas(gcf, 'scatter_str_v_ff.png');
% 
%                 % period vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, periods_vs_fills );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title(['Period (' obj.units.name ') vs. fill factors']);
%                 savefig('period_v_ff.fig');
%                 saveas(gcf, 'period_v_ff.png');
% 
%                 % offset vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, offsets_vs_fills );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Offset vs. fill factors');
%                 savefig('offsets_v_ff.fig');
%                 saveas(gcf, 'offsets_v_ff.png');
% 
%                 % k vs. fill
%                 figure;
%                 imagesc( fill_bots, fill_tops, real(k_vs_fills) );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Real k vs. fill factors');
%                 savefig('k_real_v_ff.fig');
%                 saveas(gcf, 'k_real_v_ff.png');
%                 
%                 figure;
%                 imagesc( fill_bots, fill_tops, imag(k_vs_fills) );
%                 colorbar; set( gca, 'ydir', 'normal' );
%                 xlabel('bottom fill factor'); ylabel('top fill factor');
%                 title('Imag k vs. fill factors');
%                 savefig('k_imag_v_ff.fig');
%                 saveas(gcf, 'k_imag_v_ff.png');
%                 
%             end

            % END OLD VERSION FF TOP AND BOT
            % ----------------------------------------------------------
% 
%             % ----------------------------------------
%             % NEW VERSION SWEEPING JELENAS DATAPOINTS
%             % ----------------------------------------
%             
%             fprintf('Sweeping fill factors for directivity and angle...\n');
%             
%             % set fill factors and offsets
%             fill_bots           = fliplr( 0.4:0.05:0.95 );
%             fill_top_bot_ratio  = fliplr( 0.2:0.05:1.3 );
%             fill_tops           = []; %fill_bots .* fill_top_bot_ratio;
%             offsets             = fliplr(0:0.01:0.99);
%             offsets_orig        = offsets;
%             
%             % save fills and offsets
%             obj.fill_tops           = fill_tops;
%             obj.fill_bots           = fill_bots;
%             obj.fill_top_bot_ratio  = fill_top_bot_ratio;
%             obj.offsets             = offsets;
%             
%             % initialize saving variables
%             directivities_vs_fills  = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             angles_vs_fills         = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             periods_vs_fills        = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             offsets_vs_fills        = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             scatter_str_vs_fills    = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             k_vs_fills              = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             GC_vs_fills             = cell( length( fill_bots ), length( fill_top_bot_ratio ) );      % dimensions bot fill vs. top/bot ratio
%             dir_b4_period_vs_fills  = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
%             
%             % set solver settings
%             num_modes   = 5;
%             BC          = 0;                                                % 0 = PEC
%             pml_options = [1, 200, 20, 2]; 
%             
%             tic;
%             ii = 0;
%             
%             % sweep
%             for i_ff_bot = 1:length( fill_bots )
%                 % For each bottom fill factor
%                 
%                 for i_ff_ratio = 1:length( fill_top_bot_ratio )
%                     % for each top/bottom fill factor ratio
%                     
%                     % print iteration
%                     ii = ii + 1;
%                     fprintf('Fill factor iteration %i of %i\n', ii, length( fill_top_bot_ratio ) * length( fill_bots ) );
%                     fprintf('Bottom fill factor %f of %f\n', fill_bots(i_ff_bot), fill_bots(end) );
%                     fprintf('Top/Bottom fill ratio %f of %f\n', fill_top_bot_ratio(i_ff_ratio), fill_top_bot_ratio(end) );
% 
%                     % init saving variables
%                     directivities = zeros( size(offsets) );
%                     k_vs_offset   = zeros( size(offsets) );
%                     angles        = zeros( size(offsets) );     % for debugging
% 
% 
%                     % Sweep offsets, pick offset with best directivity
%                     fprintf('Sweeping offsets...\n');
%                     for i_offset = 1:length( offsets )
% 
%                         fprintf('Iteration %i of %i\n', i_offset, length(offsets) );
% 
%                         % make grating cell
%                         fill_top = fill_top_bot_ratio(i_ff_ratio) * fill_bots(i_ff_bot);
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     guess_period, ...
%                                                     fill_top, ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     offsets(i_offset) );
%                                                 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
% 
%                         % save directivity
%                         if strcmp( obj.coupling_direction, 'up' )
%                             % coupling direction is upwards
%                             directivities( i_offset )   = GC.directivity;
%                             angles( i_offset )          = GC.max_angle_up;
%                         else
%                             % coupling direction is downwards
%                             directivities( i_offset )   = 1./( GC.directivity );
%                             angles( i_offset )          = GC.max_angle_down;
%                         end
% 
%                         % update the guessk (units rad/'units')
%                         guessk                  = GC.k;
%                         k_vs_offset( i_offset ) = GC.k;
% 
%                         toc;
% 
%                     end     % end for i_offset = ...
%                     fprintf('...done.\n');
% 
% %                         % DEBUG plot directivity vs. offset
%                     if DEBUG
%                         figure;
%                         plot( offsets, directivities, '-o' );
%                         xlabel('offsets'); ylabel('directivities');
%                         title('DEBUG directivities vs offsets');
%                         makeFigureNice();
% %                             
% %                             figure;
% %                             plot( offsets, angles, '-o' );
% %                             xlabel('offsets'); ylabel('angles');
% %                             title('DEBUG angles vs offsets for first run');
% %                             makeFigureNice();
% %                             
%                     end
% 
%                     % pick best offset
%                     [ ~, indx_best_offset ]     = max( directivities );
%                     best_offset                 = offsets( indx_best_offset );
%                     best_offset_k               = k_vs_offset( indx_best_offset );
%                     
%                     dir_b4_period_vs_fills( i_ff_bot, i_ff_ratio )      = max( directivities );
% 
%                     % DEBUG plot grating with best directivity
%                     if DEBUG
% 
%                         % make grating cell
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     guess_period, ...
%                                                     fill_tops(i_ff_top), ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     best_offset );
% 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, best_offset_k );
% 
%                         % plot field
%                         GC.plotEz_w_edges();
% 
%                     end
% 
% 
%                     % now sweep periods
%                     % only sweep larger periods. Doubtful that the period
%                     % will be smaller
%                     periods     = guess_period : obj.discretization : 1.05 * guess_period;
%                     periods     = obj.discretization * round(periods/obj.discretization);
% %                         periods_nm  = periods * obj.units.scale * 1e9;                            % convert to nm
% 
%                     % init saving variables
%                     angles          = zeros( size(periods) );
%                     k_vs_period     = zeros( size(periods) );
%                     GC_vs_period    = cell( size(periods) );
% 
%                     % sweep periods
%                     guessk = best_offset_k;
%                     fprintf('Sweeping periods...\n');
%                     for i_period = 1:length(periods)
% 
%                         fprintf('Iteration %i of %i\n', i_period, length(periods) );
% 
%                         % make grating cell
%                         fill_top = fill_top_bot_ratio(i_ff_ratio) * fill_bots(i_ff_bot);
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     periods(i_period), ...
%                                                     fill_top, ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     best_offset );
% 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
%                         % save angle
%                         if strcmp( obj.coupling_direction, 'up' )
%                             % coupling direction is upwards
%                             angles( i_period ) = GC.max_angle_up;
%                         else
%                             % coupling direction is downwards
%                             angles( i_period ) = GC.max_angle_down;
%                         end
% 
%                         % update GC list
%                         GC_vs_period{i_period} = GC;
% 
%                         % update k (units of rad/'units')
%                         k_vs_period(i_period)   = GC.k;
%                         guessk                  = GC.k;
% 
%                         toc;
% 
%                     end
%                     fprintf('...done.\n');
% 
%                     % pick best period
%                     [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles ) );
%                     best_period                     = periods( indx_best_period );
%                     best_period_k                   = k_vs_period( indx_best_period );
%                     best_GC                         = GC_vs_period{ indx_best_period };
% 
%                     
%                     
%                     % save data
%                     if strcmp( obj.coupling_direction, 'up' )
%                         % coupling direction is upwards
%                         directivities_vs_fills( i_ff_bot, i_ff_ratio )   = best_GC.directivity;
%                         angles_vs_fills( i_ff_bot, i_ff_ratio )          = best_GC.max_angle_up;
%                         scatter_str_vs_fills( i_ff_bot, i_ff_ratio )     = best_GC.alpha_up;
%                     else
%                         % coupling direction is downwards
%                         directivities_vs_fills( i_ff_bot, i_ff_ratio )   = 1./best_GC.directivity;
%                         angles_vs_fills( i_ff_bot, i_ff_ratio )          = best_GC.max_angle_down;
%                         scatter_str_vs_fills( i_ff_bot, i_ff_ratio )     = best_GC.alpha_down;
%                     end
%                     periods_vs_fills( i_ff_bot, i_ff_ratio )  = best_period;
%                     offsets_vs_fills( i_ff_bot, i_ff_ratio )  = best_offset;
%                     k_vs_fills( i_ff_bot, i_ff_ratio )        = best_GC.k;
%                     GC_vs_fills{ i_ff_bot, i_ff_ratio }       = best_GC;
%                     
%                     
%                     % update the guess parameters, period, k, offset
%                     if i_ff_ratio == 1
%                         % first iteration, save these guess parameters for
%                         % the next top level loop
%                         next_top_loop_period    = best_period;
%                         next_top_loop_k         = best_GC.k;
%                         next_top_loop_offset    = best_offset;
%                     end
%                     guessk              = best_GC.k;
%                     guess_period        = best_period;
%                     guess_offset        = best_offset;
%                     
%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                     
%                     
%                 end     % end for i_ff_ratio = ...
%                 
%                 % update the guess parameters, period, k, offset
%                 guess_period    = next_top_loop_period;
%                 guessk          = next_top_loop_k;
%                 guess_offset    = next_top_loop_offset;
%                 
%                 % update the offsets
%                 % grab previous offset index
%                 [~, indx_prev_offset] = min( abs( offsets_orig - offsets_vs_fills( i_ff_bot, 1 ) ) );
%                 % shift offsets to start at previous offset
%                 offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                 
%             end     % end for i_ff_bot = ...
%             
%             % save variables to object
%             obj.directivities_vs_fills  = directivities_vs_fills;
%             obj.angles_vs_fills         = angles_vs_fills;
%             obj.scatter_str_vs_fills    = scatter_str_vs_fills;
%             obj.periods_vs_fills        = periods_vs_fills;
%             obj.offsets_vs_fills        = offsets_vs_fills;
%             obj.k_vs_fills              = k_vs_fills;
%             obj.GC_vs_fills             = GC_vs_fills;
%             obj.dir_b4_period_vs_fills  = dir_b4_period_vs_fills;
%             

            
            % ----------------------------------------
            % NEW NEW VERSION SWEEPING JELENAS DATAPOINTS
            % SPLITTING NORMAL AND INVERTED DOMAIN
            % ----------------------------------------
            
            fprintf('Sweeping fill factors for directivity and angle...\n');
            
            % set fill factors and offsets
            fill_bots           = fliplr( 0.4:0.025:1.0 );
%             fill_top_bot_ratio  = fliplr( 0.0:0.05:1.1 );
            fill_top_bot_ratio  = fliplr( 0.0:0.025:1.2 );
%             fill_top_bot_ratio  = 1:-0.05:0.9;
%             fill_top_bot_ratio  = 0.2:-0.05:0;
%             fill_top_bot_ratio  = fliplr( 0.8:0.025:1.2 );
%             fill_bots           = 0.975;
            % on normal 1:1 line
%             fill_bots           = fliplr( 0.9:0.025:1.0 );
%             fill_top_bot_ratio  = 1;
            fill_tops           = []; %fill_bots .* fill_top_bot_ratio;
%             offsets             = fliplr(0:0.01:0.99);
%             offsets_orig        = offsets;
            guess_offset          = 0;
            
            % save fills and offsets
            obj.fill_tops           = fill_tops;
            obj.fill_bots           = fill_bots;
            obj.fill_top_bot_ratio  = fill_top_bot_ratio;
%             obj.offsets             = offsets;
            
            % split domain into inverted and normal domains
            inv_norm_thresh         = -0.79;
            fill_top_bot_ratio_inv  = fliplr( fill_top_bot_ratio( fill_top_bot_ratio < inv_norm_thresh ) );         % inverted. NOTE this array is monotonically increasing
            fill_top_bot_ratio_norm = fill_top_bot_ratio( fill_top_bot_ratio >= inv_norm_thresh );                  % normal
            
            % initialize saving variables
            % normal domain
            directivities_vs_fills_norm  = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            angles_vs_fills_norm         = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            periods_vs_fills_norm        = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            offsets_vs_fills_norm        = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio, this is offset ratio
            scatter_str_vs_fills_norm    = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            k_vs_fills_norm              = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            GC_vs_fills_norm             = cell( length( fill_bots ), length( fill_top_bot_ratio_norm ) );      % dimensions bot fill vs. top/bot ratio
            dir_b4_period_vs_fills_norm  = zeros( length( fill_bots ), length( fill_top_bot_ratio_norm ) );     % dimensions bot fill vs. top/bot ratio
            % inverted design
            directivities_vs_fills_inv  = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            angles_vs_fills_inv         = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            periods_vs_fills_inv        = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            offsets_vs_fills_inv        = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio, this is offset ratio
            scatter_str_vs_fills_inv    = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            k_vs_fills_inv              = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            GC_vs_fills_inv             = cell( length( fill_bots ), length( fill_top_bot_ratio_inv ) );      % dimensions bot fill vs. top/bot ratio
            dir_b4_period_vs_fills_inv  = zeros( length( fill_bots ), length( fill_top_bot_ratio_inv ) );     % dimensions bot fill vs. top/bot ratio
            
            
            % set solver settings
            num_modes   = 5;
            BC          = 0;                                                % 0 = PEC
            pml_options = [1, 200, 20, 2]; 
            sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);

            tic;
            ii = 0;
            
 
            
            % sweep normal
            
            % make grating cell, normal
            waveguide = obj.h_makeGratingCell( obj.convertObjToStruct(), 2*obj.discretization, 1.0, 1.0, 0.0 );
            
            % run waveguide simulation
            % sim settings
            guess_n             = 0.7 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
            guessk              = guess_n * 2*pi/obj.lambda;                                        % units rad/'units'
            num_wg_modes        = 5;
            BC                  = 0;                                                                % 0 = PEC
            pml_options_wg      = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
            % run sim
            waveguide   = waveguide.runSimulation( num_wg_modes, BC, pml_options_wg, guessk );
            
            % update guessk (units rad/'units')
            guessk = waveguide.k;
            
            % grab waveguide k
            waveguide_k = waveguide.k;  % units of rad/'units'          % * obj.units.scale * 1e9;                              % in units 1/'units'                          
            
            % DEBUG plot stuff
            waveguide.plotEz_w_edges();
            title('DEBUG waveguide, thick (normal)');
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0      = obj.background_index * ( 2*pi/obj.lambda );
            kx      = k0 * sin( (pi/180) * obj.optimal_angle );
            period  = 2*pi/(waveguide_k- kx);                                               % units of 'units'
            
            % snap period to discretization
            guess_period    = obj.discretization * round(period/obj.discretization);
            
            % ugh this is really annoying but i have to - extend the
            % waveguide's e z overlap
            [ waveguide, e_z_overlap_ext ]  = waveguide.stitch_E_field( waveguide.Phi, real(waveguide.k), round(guess_period/waveguide.domain_size(2)) );
            waveguide.E_z_for_overlap       = e_z_overlap_ext;
            
            fprintf('...done\n\n');
            
            % now run the sweep
            
            % only run if normal domain exists
            if length(fill_top_bot_ratio_norm) > 0
                
                % initially start with waveguide GC
                guess_GC = waveguide;
                
                % first fill in the right side of the domain
                for i_ff_bot = 1:length( fill_bots )

                    % print iteration
                    fprintf('Right side of domain, iteration %i of %i\n', i_ff_bot, length(fill_bots) );

                    % Optimize period and offset
                    fill_top = fill_top_bot_ratio_norm(1) * fill_bots(i_ff_bot);
                    if fill_top < 1
                        % Only run optimization if theres a perturbation 
                        
                        [ obj, best_period, best_offset, best_directivity, ...
                          best_angle, best_scatter_str, best_GC, ...
                          best_k, dir_b4_period_vs_fill ] = obj.optimizePeriodOffset( guess_offset, ...
                                                                                      fill_top, ...
                                                                                      fill_bots(i_ff_bot), ...
                                                                                      guess_period,...
                                                                                      guessk, ...
                                                                                      sim_opts, ...
                                                                                      guess_GC );

                        % save data
                        if strcmp( obj.coupling_direction, 'up' )
                            % coupling direction is upwards
                            directivities_vs_fills_norm( i_ff_bot, 1 )   = best_GC.directivity;
                            angles_vs_fills_norm( i_ff_bot, 1 )          = best_GC.max_angle_up;
                            scatter_str_vs_fills_norm( i_ff_bot, 1 )     = best_GC.alpha_up;
                        else
                            % coupling direction is downwards
                            directivities_vs_fills_norm( i_ff_bot, 1 )   = 1./best_GC.directivity;
                            angles_vs_fills_norm( i_ff_bot, 1 )          = best_GC.max_angle_down;
                            scatter_str_vs_fills_norm( i_ff_bot, 1 )     = best_GC.alpha_down;
                        end
                        periods_vs_fills_norm( i_ff_bot, 1 )          = best_period;
                        offsets_vs_fills_norm( i_ff_bot, 1 )          = best_offset/best_period;
                        k_vs_fills_norm( i_ff_bot, 1 )                = best_GC.k;
                        GC_vs_fills_norm{ i_ff_bot, 1 }               = best_GC;
                        dir_b4_period_vs_fills_norm( i_ff_bot, 1 )    = dir_b4_period_vs_fill;


                        % update the guess parameters, period, k, offset
                        guessk              = best_GC.k;
                        guess_period        = best_period;
                        guess_GC            = best_GC;
                        guess_offset        = best_offset;

%                         % update the offsets
%                         % grab previous offset index
%                         [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                         % shift offsets to start at previous offset
%                         offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
                        
                    else
                        % we're in a waveguide, there's no reason to run
                        % the optimization (and actually the period sweep
                        % bugs out when the fill = 100%)
                        
                        % save dummy 
                        directivities_vs_fills_norm( i_ff_bot, 1 )   = 1;
                        angles_vs_fills_norm( i_ff_bot, 1 )          = 0;
                        scatter_str_vs_fills_norm( i_ff_bot, 1 )     = 0;
                        periods_vs_fills_norm( i_ff_bot, 1 )          = guess_period;
                        offsets_vs_fills_norm( i_ff_bot, 1 )          = 0;
                        k_vs_fills_norm( i_ff_bot, 1 )                = guessk;
                        GC_vs_fills_norm{ i_ff_bot, 1 }               = waveguide;
                        dir_b4_period_vs_fills_norm( i_ff_bot, 1 )    = 1;
                        
                    end

                    toc;

                end     % end initial normal domain sweep


                % for parallel processing, grab these variables before entering
                % parfor
                offsets_vs_fills_norm_1 = offsets_vs_fills_norm(:,1);
                periods_vs_fills_norm_1 = periods_vs_fills_norm(:,1);
                k_vs_fills_norm_1       = k_vs_fills_norm(:,1);
                GC_vs_fills_norm_1      = GC_vs_fills_norm(:,1);
                % calc number of loops, also necessary apparently for parfor
                n_fill_bots                 = length(fill_bots);
                n_fill_top_bot_ratio_norm   = length(fill_top_bot_ratio_norm);
                
                % start a parallel pool session
                my_cluster = parcluster('local');                       % cores on compute node are "local"
                if getenv('ENVIRONMENT')                                % true if this is a batch job
                    my_cluster.JobStorageLocation = getenv('TMPDIR');    % points to TMPDIR
                end

                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if ~isempty(poolobj)
                    % shut down previously made parallel pool
                    delete(gcp('nocreate'));
                end
                parpool(my_cluster, my_cluster.NumWorkers);


                % create a parforprogmon object
%                 ppm = ParforProgMon('Progress on normal parfor: ', n_fill_bots );

                % now fill in the rest of the domain
                parfor i_ff_bot = 1:n_fill_bots
                    % For each bottom fill factor

    %                 % skip this loop if we aren't running the normal domain
    %                 if length( fill_top_bot_ratio_norm ) == 0
    %                     break;
    %                 end
    
                    fprintf('Normal parfor iteration %i of %i\n', i_ff_bot, n_fill_bots);

%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - offsets_vs_fills_norm_1( i_ff_bot ) ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );

                    % grab starting guess period and k
                    guess_period    = periods_vs_fills_norm_1( i_ff_bot );
                    guessk          = k_vs_fills_norm_1( i_ff_bot );
                    guess_GC        = GC_vs_fills_norm_1{ i_ff_bot };
                    guess_offset    = offsets_vs_fills_norm_1( i_ff_bot ) * guess_period;

                    % grab bottom fill
                    fill_bot = fill_bots( i_ff_bot );


                    for i_ff_ratio_norm = 2:n_fill_top_bot_ratio_norm
                        % for each top/bottom fill factor ratio

                         fprintf('Fill factor ratio %i of %i, Normal parfor iteration %i of %i\n', i_ff_ratio_norm, n_fill_top_bot_ratio_norm, i_ff_bot, n_fill_bots);
                        
                        % Optimize period and offset
                        fill_top = fill_top_bot_ratio_norm(i_ff_ratio_norm) * fill_bot;
                        if fill_top < 1
                            % Only run optimization if theres a perturbation 
                            [ tempobj, best_period, best_offset, best_directivity, ...
                              best_angle, best_scatter_str, best_GC, ...
                              best_k, dir_b4_period_vs_fill ] = obj.optimizePeriodOffset( guess_offset, ...
                                                                                          fill_top, ...
                                                                                          fill_bot, ...
                                                                                          guess_period,...
                                                                                          guessk, ...
                                                                                          sim_opts, ...
                                                                                          guess_GC );

                            % save data
                            if strcmp( obj.coupling_direction, 'up' )
        %                         coupling direction is upwards
                                directivities_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )   = best_GC.directivity;
                                angles_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = best_GC.max_angle_up;
                                scatter_str_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )     = best_GC.alpha_up;
                            else
                                % coupling direction is downwards
                                directivities_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )   = 1./best_GC.directivity;
                                angles_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = best_GC.max_angle_down;
                                scatter_str_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )     = best_GC.alpha_down;
                            end
                            periods_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = best_period;
                            offsets_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = best_offset/best_period;
                            k_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )                = best_GC.k;
                            GC_vs_fills_norm{ i_ff_bot, i_ff_ratio_norm }               = best_GC;
                            dir_b4_period_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )    = dir_b4_period_vs_fill;

                            % update the guess parameters, period, k, offset
                            guessk              = best_GC.k;
                            guess_period        = best_period;
                            guess_GC            = best_GC;
                            guess_offset        = best_offset;

%                             % update the offsets
%                             % grab previous offset index
%                             [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                             % shift offsets to start at previous offset
%                             offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
                            
                        else
                            % we're in a waveguide, there's no reason to run
                            % the optimization (and actually the period sweep
                            % bugs out when the fill = 100%)

                            % save dummy 
                            directivities_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )    = 1;
                            angles_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )           = 0;
                            scatter_str_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )      = 0;
                            periods_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = guess_period;
                            offsets_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )          = 0;
                            k_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )                = guessk;
                            GC_vs_fills_norm{ i_ff_bot, i_ff_ratio_norm }               = waveguide;
                            dir_b4_period_vs_fills_norm( i_ff_bot, i_ff_ratio_norm )    = 1;

                        end


                    end     % end for i_ff_ratio = ...

                    % update progress monitor
%                     ppm.increment();

                end     % end for i_ff_bot = ...
                
            end     % end if length(fill_top_bot_ratio_norm) > 0
            
            
            
            % sweep inverted
            
            % make grating cell, inverted
            waveguide = obj.h_makeGratingCell( obj.convertObjToStruct(), obj.discretization, 0.0, 1.0, 0.0 );
            
            % run waveguide simulation
            % sim settings
            guess_n         = 0.7 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
            guessk          = guess_n * 2*pi/obj.lambda;                                        % units rad/'units'
            num_wg_modes    = 5;
            BC              = 0;                                                                % 0 = PEC
            pml_options_wg  = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
            % run sim
            waveguide   = waveguide.runSimulation( num_wg_modes, BC, pml_options_wg, guessk );
            
            % update guessk (units rad/'units')
            guessk = waveguide.k;
            
            % grab waveguide k
            waveguide_k = waveguide.k;  % units of rad/'units'          % * obj.units.scale * 1e9;                              % in units 1/'units'                          
            
            % DEBUG plot stuff
            waveguide.plotEz_w_edges();
            title('DEBUG waveguide, thin (inverted)');
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0      = obj.background_index * ( 2*pi/obj.lambda );
            kx      = k0 * sin( (pi/180) * obj.optimal_angle );
            period  = 2*pi/(waveguide_k - kx);                                               % units of 'units'
            
            % snap period to discretization
            guess_period    = obj.discretization * round(period/obj.discretization);
    
            % ugh this is really annoying but i have to - extend the
            % waveguide's e z overlap
            [ waveguide, e_z_overlap_ext ]  = waveguide.stitch_E_field( waveguide.Phi, real(waveguide.k), round(guess_period/waveguide.domain_size(2)) );
            waveguide.E_z_for_overlap       = e_z_overlap_ext;
            
            fprintf('...done\n\n');
            
            % reset offsets 
            guess_offset = 0;
            
            % now run the sweep
            
            % only run if invert domain exists
            if length(fill_top_bot_ratio_inv) > 0
            
                % initially start with waveguide GC
                guess_GC = waveguide;
                
                % first fill in the left side of the domain
                for i_ff_bot = 1:length( fill_bots )

                    % print iteration
                    fprintf('Left side of domain, iteration %i of %i\n', i_ff_bot, length(fill_bots) );

                    % Optimize period and offset

                    fill_top = fill_top_bot_ratio_inv(1) * fill_bots(i_ff_bot);
                    [ obj, best_period, best_offset, best_directivity, ...
                      best_angle, best_scatter_str, best_GC, ...
                      best_k, dir_b4_period_vs_fill ] = obj.optimizePeriodOffset( guess_offset, ...
                                                                                  fill_top, ...
                                                                                  fill_bots(i_ff_bot), ...
                                                                                  guess_period,...
                                                                                  guessk, ...
                                                                                  sim_opts, ...
                                                                                  guess_GC );


                    % save data
                    if strcmp( obj.coupling_direction, 'up' )
                        % coupling direction is upwards
                        directivities_vs_fills_inv( i_ff_bot, 1 )   = best_GC.directivity;
                        angles_vs_fills_inv( i_ff_bot, 1 )          = best_GC.max_angle_up;
                        scatter_str_vs_fills_inv( i_ff_bot, 1 )     = best_GC.alpha_up;
                    else
                        % coupling direction is downwards
                        directivities_vs_fills_inv( i_ff_bot, 1 )   = 1./best_GC.directivity;
                        angles_vs_fills_inv( i_ff_bot, 1 )          = best_GC.max_angle_down;
                        scatter_str_vs_fills_inv( i_ff_bot, 1 )     = best_GC.alpha_down;
                    end
                    periods_vs_fills_inv( i_ff_bot, 1 )          = best_period;
                    offsets_vs_fills_inv( i_ff_bot, 1 )          = best_offset/best_period;
                    k_vs_fills_inv( i_ff_bot, 1 )                = best_GC.k;
                    GC_vs_fills_inv{ i_ff_bot, 1 }               = best_GC;
                    dir_b4_period_vs_fills_inv( i_ff_bot, 1 )    = dir_b4_period_vs_fill;


                    % update the guess parameters, period, k, offset
                    guessk              = best_GC.k;
                    guess_period        = best_period;
                    guess_GC            = best_GC;
                    guess_offset        = best_offset;

%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );

                    toc;

                end     % end initial normal domain sweep

                % for parallel processing, grab these variables before entering
                % parfor
                offsets_vs_fills_inv_1 = offsets_vs_fills_inv(:,1);
                periods_vs_fills_inv_1 = periods_vs_fills_inv(:,1);
                k_vs_fills_inv_1       = k_vs_fills_inv(:,1);
                GC_vs_fills_inv_1      = GC_vs_fills_inv(:,1);
                % calc number of loops, also necessary apparently for parfor
                n_fill_bots                 = length(fill_bots);
                n_fill_top_bot_ratio_inv    = length(fill_top_bot_ratio_inv);

                
                % start a parallel pool session
                my_cluster = parcluster('local');                       % cores on compute node are "local"
                if getenv('ENVIRONMENT')                                % true if this is a batch job
                    my_cluster.JobStorageLocation = getenv('TMPDIR');    % points to TMPDIR
                end

                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if ~isempty(poolobj)
                    % shut down previously made parallel pool
                    delete(gcp('nocreate'));
                end
                parpool(my_cluster, my_cluster.NumWorkers);

                % create a parforprogmon object
%                 ppm = ParforProgMon('Progress on inverted parfor: ', n_fill_bots );

                % now fill in the rest of the domain
                parfor i_ff_bot = 1:n_fill_bots
                    % For each bottom fill factor

                    fprintf('Invert parfor iteration %i of %i\n', i_ff_bot, n_fill_bots);
                    
%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - offsets_vs_fills_inv_1( i_ff_bot ) ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );

                    % grab starting guess period and k
                    guess_period    = periods_vs_fills_inv_1( i_ff_bot );
                    guessk          = k_vs_fills_inv_1( i_ff_bot );
                    guess_GC        = GC_vs_fills_inv_1{ i_ff_bot };
                    guess_offset    = offsets_vs_fills_inv_1( i_ff_bot );

                    % grab bottom fill
                    fill_bot = fill_bots( i_ff_bot );


                    for i_ff_ratio_inv = 2:n_fill_top_bot_ratio_inv
                        % for each top/bottom fill factor ratio

                        % Optimize period and offset
                        fill_top = fill_top_bot_ratio_inv(i_ff_ratio_inv) * fill_bot;
                        [ tempobj, best_period, best_offset, best_directivity, ...
                          best_angle, best_scatter_str, best_GC, ...
                          best_k, dir_b4_period_vs_fill ] = obj.optimizePeriodOffset( guess_offset, ...
                                                                                      fill_top, ...
                                                                                      fill_bot, ...
                                                                                      guess_period,...
                                                                                      guessk, ...
                                                                                      sim_opts, ...
                                                                                      guess_GC );


                        % save data
                        if strcmp( obj.coupling_direction, 'up' )
                            % coupling direction is upwards
                            directivities_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )   = best_GC.directivity;
                            angles_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_GC.max_angle_up;
                            scatter_str_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )     = best_GC.alpha_up;
                        else
                            % coupling direction is downwards
                            directivities_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )   = 1./best_GC.directivity;
                            angles_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_GC.max_angle_down;
                            scatter_str_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )     = best_GC.alpha_down;
                        end
                        periods_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_period;
                        offsets_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_offset/best_period;
                        k_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )                = best_GC.k;
                        GC_vs_fills_inv{ i_ff_bot, i_ff_ratio_inv }               = best_GC;
                        dir_b4_period_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )    = dir_b4_period_vs_fill;


                        % update the guess parameters, period, k, offset
                        guessk              = best_GC.k;
                        guess_period        = best_period;
                        guess_GC            = best_GC;
                        guess_offset        = best_offset;

%                         % update the offsets
%                         % grab previous offset index
%                         [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                         % shift offsets to start at previous offset
%                         offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );


                    end     % end for i_ff_ratio = ...

                    % update progress monitor
%                     ppm.increment();

                end     % end for i_ff_bot = ...
                
            end     % end if length(fill_top_bot_ratio_inv) > 0
            
            % OLD LOOP
%             for i_ff_bot = 1:length( fill_bots )
%                 % For each bottom fill factor
%                 
% %                 if length( fill_top_bot_ratio_inv ) == 0
% %                     break;
% %                 end
%                 
%                 for i_ff_ratio_inv = 1:length( fill_top_bot_ratio_inv )
%                     % for each top/bottom fill factor ratio
%                     
%                     % print iteration
%                     ii = ii + 1;
%                     fprintf('Fill factor iteration %i of %i\n', ii, length( fill_top_bot_ratio ) * length( fill_bots ) );
%                     fprintf('Bottom fill factor %f of %f\n', fill_bots(i_ff_bot), fill_bots(end) );
%                     fprintf('Top/Bottom fill ratio %f of %f\n', fill_top_bot_ratio_inv(i_ff_ratio_inv), fill_top_bot_ratio_inv(end) );
% 
%                     % init saving variables
%                     directivities = zeros( size(offsets) );
%                     k_vs_offset   = zeros( size(offsets) );
%                     angles        = zeros( size(offsets) );     % for debugging
% %                     GC_vs_offset  = cell( size(offsets) );
% 
% 
%                     % Sweep offsets, pick offset with best directivity
%                     fprintf('Sweeping offsets...\n');
%                     for i_offset = 1:length( offsets )
% 
%                         fprintf('Iteration %i of %i\n', i_offset, length(offsets) );
% 
%                         % make grating cell
%                         fill_top = fill_top_bot_ratio_inv(i_ff_ratio_inv) * fill_bots(i_ff_bot);
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     guess_period, ...
%                                                     fill_top, ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     offsets(i_offset) );
%                                                 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
% 
%                         % save directivity
%                         if strcmp( obj.coupling_direction, 'up' )
%                             % coupling direction is upwards
%                             directivities( i_offset )   = GC.directivity;
%                             angles( i_offset )          = GC.max_angle_up;
%                         else
%                             % coupling direction is downwards
%                             directivities( i_offset )   = 1./( GC.directivity );
%                             angles( i_offset )          = GC.max_angle_down;
%                         end
% 
%                         % update the guessk (units rad/'units')
%                         guessk                  = GC.k;
%                         k_vs_offset( i_offset ) = GC.k;
% %                         GC_vs_offset{ i_offset } = GC;
% 
%                         toc;
% 
%                     end     % end for i_offset = ...
%                     fprintf('...done.\n');
% 
% %                         % DEBUG plot directivity vs. offset
%                     if DEBUG
%                         figure;
%                         plot( offsets, directivities, '-o' );
%                         xlabel('offsets'); ylabel('directivities');
%                         title('DEBUG directivities vs offsets');
%                         makeFigureNice();
% %                             
% %                             figure;
% %                             plot( offsets, angles, '-o' );
% %                             xlabel('offsets'); ylabel('angles');
% %                             title('DEBUG angles vs offsets for first run');
% %                             makeFigureNice();
% %                             
%                     end
% 
%                     % pick best offset
%                     [ ~, indx_best_offset ]     = max( directivities );
%                     best_offset                 = offsets( indx_best_offset );
%                     best_offset_k               = k_vs_offset( indx_best_offset );
%                     best_offset_angle           = angles( indx_best_offset );
% %                     best_offset_GC              = GC_vs_offset{ indx_best_offset };
%                     
%                     dir_b4_period_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )      = max( directivities );
% 
%                     % DEBUG plot grating with best directivity
%                     if DEBUG
% 
%                         % make grating cell
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     guess_period, ...
%                                                     fill_tops(i_ff_top), ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     best_offset );
% 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, best_offset_k );
% 
%                         % plot field
%                         GC.plotEz_w_edges();
% 
%                     end
% 
% 
%                     % now sweep periods
%                     
%                     % decide whether to sweep larger or smaller periods
%                     % based on the angle
%                     if best_offset_angle > obj.optimal_angle
%                         % only sweep smaller periods
%                         periods     = guess_period : -obj.discretization : 0.90 * guess_period;
%                         periods     = obj.discretization * round(periods/obj.discretization);
%                     else
%                         % only sweep larger periods
%                         periods     = guess_period : obj.discretization : 1.05 * guess_period;
%                         periods     = obj.discretization * round(periods/obj.discretization);
%                     end
% 
%                     % init saving variables
%                     angles          = zeros( size(periods) );
%                     k_vs_period     = zeros( size(periods) );
%                     GC_vs_period    = cell( size(periods) );
% 
%                     % sweep periods
%                     guessk = best_offset_k;
%                     fprintf('Sweeping periods...\n');
%                     for i_period = 1:length(periods)
% 
%                         fprintf('Iteration %i of %i\n', i_period, length(periods) );
% 
%                         % make grating cell
%                         fill_top = fill_top_bot_ratio_inv(i_ff_ratio_inv) * fill_bots(i_ff_bot);
%                         GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                                     periods(i_period), ...
%                                                     fill_top, ...
%                                                     fill_bots(i_ff_bot), ...
%                                                     best_offset );
% 
%                         % run sim
%                         GC = GC.runSimulation( num_modes, BC, pml_options, guessk );
% 
%                         % save angle
%                         if strcmp( obj.coupling_direction, 'up' )
%                             % coupling direction is upwards
%                             angles( i_period ) = GC.max_angle_up;
%                         else
%                             % coupling direction is downwards
%                             angles( i_period ) = GC.max_angle_down;
%                         end
% 
%                         % update GC list
%                         GC_vs_period{i_period} = GC;
% 
%                         % update k (units of rad/'units')
%                         k_vs_period(i_period)   = GC.k;
%                         guessk                  = GC.k;
% 
%                         toc;
% 
%                     end
%                     fprintf('...done.\n');
% 
%                     % pick best period
%                     [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles ) );
%                     best_period                     = periods( indx_best_period );
%                     best_period_k                   = k_vs_period( indx_best_period );
%                     best_GC                         = GC_vs_period{ indx_best_period };
% 
%                     
%                     
%                     % save data
%                     if strcmp( obj.coupling_direction, 'up' )
%                         % coupling direction is upwards
%                         directivities_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )   = best_GC.directivity;
%                         angles_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_GC.max_angle_up;
%                         scatter_str_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )     = best_GC.alpha_up;
%                     else
%                         % coupling direction is downwards
%                         directivities_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )   = 1./best_GC.directivity;
%                         angles_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )          = best_GC.max_angle_down;
%                         scatter_str_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )     = best_GC.alpha_down;
%                     end
%                     periods_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )  = best_period;
%                     offsets_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )  = best_offset;
%                     k_vs_fills_inv( i_ff_bot, i_ff_ratio_inv )        = best_GC.k;
%                     GC_vs_fills_inv{ i_ff_bot, i_ff_ratio_inv }       = best_GC;
%                     
%                     
%                     % update the guess parameters, period, k, offset
%                     if i_ff_ratio_inv == 1
%                         % first iteration, save these guess parameters for
%                         % the next top level loop
%                         next_top_loop_period    = best_period;
%                         next_top_loop_k         = best_GC.k;
%                         next_top_loop_offset    = best_offset;
%                     end
%                     guessk              = best_GC.k;
%                     guess_period        = best_period;
%                     guess_offset        = best_offset;
%                     
%                     % update the offsets
%                     % grab previous offset index
%                     [~, indx_prev_offset] = min( abs( offsets_orig - best_offset ) );
%                     % shift offsets to start at previous offset
%                     offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                     
%                     
%                 end     % end for i_ff_ratio = ...
%                 
%                 % update the guess parameters, period, k, offset
%                 guess_period    = next_top_loop_period;
%                 guessk          = next_top_loop_k;
%                 guess_offset    = next_top_loop_offset;
%                 
%                 % update the offsets
%                 % grab previous offset index
%                 [~, indx_prev_offset] = min( abs( offsets_orig - offsets_vs_fills_inv( i_ff_bot, 1 ) ) );
%                 % shift offsets to start at previous offset
%                 offsets = circshift( offsets_orig, -( indx_prev_offset - 1 ) );
%                 
%             end     % end for i_ff_bot = ...
            
            
            
            % save variables to object
            obj.directivities_vs_fills  = [ directivities_vs_fills_norm, fliplr( directivities_vs_fills_inv ) ];
            obj.angles_vs_fills         = [ angles_vs_fills_norm, fliplr( angles_vs_fills_inv ) ];
            obj.scatter_str_vs_fills    = [ scatter_str_vs_fills_norm, fliplr( scatter_str_vs_fills_inv ) ];
            obj.periods_vs_fills        = [ periods_vs_fills_norm, fliplr( periods_vs_fills_inv ) ];
            obj.offsets_vs_fills        = [ offsets_vs_fills_norm, fliplr( offsets_vs_fills_inv ) ];
            obj.k_vs_fills              = [ k_vs_fills_norm, fliplr( k_vs_fills_inv ) ];
            obj.GC_vs_fills             = [ GC_vs_fills_norm, fliplr( GC_vs_fills_inv ) ];
            obj.dir_b4_period_vs_fills  = [ dir_b4_period_vs_fills_norm, fliplr( dir_b4_period_vs_fills_inv ) ];


            % generate final designs:
%             obj = obj.generateFinalDesignGaussian( MFD );
%             obj = obj.runFinalDesignEME( MFD );

%             % run final design in EME
%             obj = obj.runFinalDesignEME(MFD);
            
        end     % end synthesizeGaussianGrating()
        
        
        function [  obj, best_period, best_offset, best_directivity, best_angle, best_scatter_str, ...
                    best_GC, best_k, dir_b4_period_vs_fill ] ...
                    = optimizePeriodOffset(obj, guess_offset, fill_top, fill_bot, guess_period, guessk, sim_opts, guess_gc )
            % optimizes period and offset for best angle/directivity
            %
            % Inputs:
            %   offsets - OLD BUT STILL IN USE
            %   guess_offset - NEW BUT NOT IMPLEMENTED
            %       type: scalar, double
            %       desc: guess offset position to start from, in units
            %             'units'
            %   fill_top
            %       type: scalar, double
            %       desc: ratio of top waveguide to period
            %   fill_bot
            %       type: scalar, double
            %       desc: ratio of bottom waveguide to period
            %   guess_period
            %       type: scalar, double
            %       desc: starting period to sweep from, in units 'units'
            %   guessk
            %       type: scalar, double
            %       desc: starting guess k
            %   sim_opts
            %       type: struct
            %       desc: structure with these fields:  
            %               num_modes
            %               BC
            %               pml_options
            %   guess_gc
            %       type: grating coupler object
            %       desc: initial grating coupler object to start with (for mode
            %       overlapping)
            %       
            %
            % Outputs:
            %   best_period
            %       type:
            %       desc: 
            %   best_offset
            %       type:
            %       desc: absolute value of offset, in 'units'
            %   best_offset_ratio
            %       type:
            %       desc:
            %   best_directivity
            %       type:
            %       desc:
            %   best_angle
            %       type:
            %       desc:
            %   best_scatter_str
            %       type:
            %       desc:
            %   best_GC
            %       type:
            %       desc:
            %   best_k
            %       type:
            %       desc:
            %   dir_b4_period_vs_fill
            %       type:
            %       desc:
            %       mostly for debugging purposes
            
%             % enable/disable debug mode
%             DEBUG = false;
            
            % generate vector of absolute offsets
            offsets         = guess_offset : obj.discretization : guess_offset + guess_period;
            offsets         = mod( offsets, guess_period );
            offset_ratios   = offsets / guess_period;

            % init saving variables vs. offset
            directivities = zeros( size(offset_ratios) );
            k_vs_offset   = zeros( size(offset_ratios) );
            angles        = zeros( size(offset_ratios) );
            GC_vs_offset  = {};

            % grab mode to overlap with
            OPTS = struct( 'mode_to_overlap', guess_gc.E_z_for_overlap );

            % Sweep offsets, pick offset with best directivity
%             fprintf('Sweeping offsets...\n');
            for i_offset = 1:length( offset_ratios )

%                 fprintf('Iteration %i of %i\n', i_offset, length(offsets) );

                % make grating cell
%                 fill_top = fill_top_bot_ratio_norm(i_ff_ratio_norm) * fill_bots(i_ff_bot);
                GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                            guess_period, ...
                                            fill_top, ...
                                            fill_bot, ...
                                            offset_ratios(i_offset) );

                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, guessk, OPTS );

                % save directivity
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    directivities( i_offset )   = GC.directivity;
                    angles( i_offset )          = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    directivities( i_offset )   = 1./( GC.directivity );
                    angles( i_offset )          = GC.max_angle_down;
                end

                % update the guessk (units rad/'units')
                guessk                  = GC.k;
                k_vs_offset( i_offset ) = GC.k;
                
                % update the mode overlap field
                GC_vs_offset{i_offset}  = GC;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;

%                 toc;

            end     % end for i_offset = ...
%             fprintf('...done.\n');

%                         % DEBUG plot directivity vs. offset
%             if DEBUG
%                 figure;
%                 plot( offsets, directivities, '-o' );
%                 xlabel('offsets'); ylabel('directivities');
%                 title('DEBUG directivities vs offsets');
%                 makeFigureNice();
%                             
%                             figure;
%                             plot( offsets, angles, '-o' );
%                             xlabel('offsets'); ylabel('angles');
%                             title('DEBUG angles vs offsets for first run');
%                             makeFigureNice();
%                             
%             end

            % pick best offset
            [ ~, indx_best_offset ]     = max( directivities );
%             best_offset                 = offsets( indx_best_offset );
            best_offset_ratio           = offset_ratios( indx_best_offset );
            best_offset_k               = k_vs_offset( indx_best_offset );
            best_offset_angle           = angles( indx_best_offset );
            
            % update mode overlap
            OPTS.mode_to_overlap    = GC_vs_offset{indx_best_offset}.E_z_for_overlap;

            % here's an output variable
            dir_b4_period_vs_fill = max( directivities );

            % DEBUG plot grating with best directivity
%             if DEBUG
% 
%                 % make grating cell
%                 GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
%                                             guess_period, ...
%                                             fill_tops(i_ff_top), ...
%                                             fill_bots(i_ff_bot), ...
%                                             best_offset );
% 
%                 % run sim
%                 GC = GC.runSimulation( num_modes, BC, pml_options, best_offset_k );
% 
%                 % plot field
%                 GC.plotEz_w_edges();
% 
%             end


%             % now sweep periods
%             % only sweep larger periods. Doubtful that the period
%             % will be smaller
%             periods     = guess_period : obj.discretization : 1.05 * guess_period;
%             periods     = obj.discretization * round(periods/obj.discretization);
% %                         periods_nm  = periods * obj.units.scale * 1e9;                            % convert to nm

            % now sweep periods
            % decide whether to sweep larger or smaller periods
            % based on the angle
            if best_offset_angle > obj.optimal_angle
                % only sweep smaller periods
                decrease_periods = true;
            else
                % only sweep larger periods
                decrease_periods = false;
            end

            % init saving variables
            angles_vs_period    = []; %= zeros( size(periods) );
            k_vs_period         = []; %zeros( size(periods) );
            GC_vs_period        = {}; % cell( size(periods) );
            periods             = [];

            % sweep periods
            guessk = best_offset_k;
            period = guess_period;
%             fprintf('Sweeping periods...\n');
            
            % set while loop exit flag
            angle_err_sign_flip = false;
            
%             for i_period = 1:length(periods)
            i_period    = 0;
            while ~angle_err_sign_flip

                i_period = i_period + 1;
%                 fprintf('Iteration %i\n', i_period );

                % make grating cell
%                 fill_top = fill_top_bot_ratio_norm(i_ff_ratio_norm) * fill_bots(i_ff_bot);
                GC = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                            period, ...
                                            fill_top, ...
                                            fill_bot, ...
                                            best_offset_ratio );

                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, guessk, OPTS );

                % save angle
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    angles_vs_period( i_period ) = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    angles_vs_period( i_period ) = GC.max_angle_down;
                end

                % update for next iteration
                periods(i_period)       = period;
                GC_vs_period{i_period}  = GC;
                k_vs_period(i_period)   = GC.k;
                guessk                  = GC.k;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;
                
                % update period
                if decrease_periods == true
                    
                    % check for exit condition (error of angle switches
                    % sign)
                    if angles_vs_period(i_period) < obj.optimal_angle
                        angle_err_sign_flip = true;
                    else
                        % decrease the period
                        period = period - obj.discretization;
                    end
                    
                else
                    
                    % check for exit condition (error of angle switches
                    % sign)
                    if angles_vs_period(i_period) > obj.optimal_angle
                        angle_err_sign_flip = true;
                    else
                        % increase the period
                        period = period + obj.discretization;
                    end
                    
                end     % end updating period if else
                
%                 toc;

            end     % end period sweep
%             fprintf('...done.\n');

            % pick best period
            [angle_error, indx_best_period] = min( abs( obj.optimal_angle - angles_vs_period ) );
            best_period_k                   = k_vs_period( indx_best_period );
            
 
            % return data
            % best offset is already set
            best_GC = GC_vs_period{ indx_best_period };
            
            if strcmp( obj.coupling_direction, 'up' )
                % coupling direction is upwards
                best_directivity    = best_GC.directivity;
                best_angle          = best_GC.max_angle_up;
                best_scatter_str    = best_GC.alpha_up;
            else
                % coupling direction is downwards
                best_directivity    = 1./best_GC.directivity;
                best_angle          = best_GC.max_angle_down;
                best_scatter_str    = best_GC.alpha_down;
            end
            
            best_period = periods( indx_best_period );
            best_k      = best_GC.k;
            best_offset = best_offset_ratio * best_period;
                    
        end     % end function optimizePeriodOffset()
        
        
        function obj = generateFinalDesignGaussian(obj, MFD, input_wg_type)
            % function for generating the final synthesized design
            % parameters
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter, in units 'units'
            %   input_wg_type
            %       type: string
            %       desc: 'bottom' for cSi only or 'full' for both cSi and
            %             pSi
            %
            % Sets these fields:
            %   obj.dir_synth                   = [];
            %   obj.bot_fill_synth              = [];
            %   obj.top_bot_fill_ratio_synth    = [];
            %   obj.period_synth                = [];
            %   obj.offset_synth                = [];
            %   obj.angles_synth                = [];
            %   obj.scatter_str_synth           = [];
            %   obj.k_synth                     = [];
            %   obj.GC_synth                    = {};
            %   obj.des_scatter_norm            = [];
            

            % save input waveguide type
            obj.input_wg_type = input_wg_type;
            
            % generate x coordinates for the gaussian mode
            % must be large enough to fit mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiberModeGaussian(    w0, zvec, xvec,...
                                                    obj.optimal_angle, d0, obj.background_index );
            
            % calculate desired scattering strength vs. x
            integral_u      = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
            alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );  % in units 1/m
%             alpha_des(end)  = 0;                                                            % for stability
            alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
            
            
            % DEBUG plot alpha desired
            figure;
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice();
%             

            % meshgrid the fills
            [ topbot_ratio_mesh, bot_fills_mesh ] = meshgrid( obj.fill_top_bot_ratio, obj.fill_bots );
            
            % set invert/normal threshold
            invert_normal_top_bot_ratio_thresh = 0.6;

            if strcmp( input_wg_type, 'bottom' ) == true
                % Inverted design
            
                % first narrow down the space to take the datapoints with the
                % maximum directivity per bottom fill factor (so for each row
                % on my design space)

                indx_bottom_space   = obj.fill_top_bot_ratio < invert_normal_top_bot_ratio_thresh;

                % grab variables (remember dimensions are bot fill x top bot
                % ratio)
                topbot_ratio_vs_fills_inv   = topbot_ratio_mesh( :, indx_bottom_space );
                bot_fills_vs_fills_inv      = bot_fills_mesh( :, indx_bottom_space );
                directivities_vs_fills_inv  = obj.directivities_vs_fills( :, indx_bottom_space );
                angles_vs_fills_inv         = obj.angles_vs_fills( :, indx_bottom_space );      
                periods_vs_fills_inv        = obj.periods_vs_fills( :, indx_bottom_space );
                offsets_vs_fills_inv        = obj.offsets_vs_fills( :, indx_bottom_space );
                scatter_str_vs_fills_inv    = obj.scatter_str_vs_fills( :, indx_bottom_space );
                k_vs_fills_inv              = obj.k_vs_fills( :, indx_bottom_space );
                GC_vs_fills_inv             = obj.GC_vs_fills( :, indx_bottom_space );
                
                % remove datapoints where the angle deviates beyond some
                % angle, by artifically setting the directivity to be very
                % low
                directivities_vs_fills_inv( abs( angles_vs_fills_inv - obj.optimal_angle ) > 1 ) = 1e-6;


                % for each top/bottom ratio, pick bottom fill with highest
                % directivity
                [ highest_dir_per_ratio, indx_highest_dir_per_ratio ] = max( directivities_vs_fills_inv, [], 1 );
                % gotta use linear indexing
                indxs                   = sub2ind( size(angles_vs_fills_inv), indx_highest_dir_per_ratio, 1:size(angles_vs_fills_inv,2)  );
                angles_high_dir         = angles_vs_fills_inv( indxs );
                periods_high_dir        = periods_vs_fills_inv( indxs );
                offsets_high_dir        = offsets_vs_fills_inv( indxs );
                scatter_strs_high_dir   = scatter_str_vs_fills_inv( indxs );
                k_high_dir              = k_vs_fills_inv( indxs );
                topbot_ratio_high_dir   = topbot_ratio_vs_fills_inv( indxs );
                bot_fills_high_dir      = bot_fills_vs_fills_inv( indxs );
                GC_high_dir             = GC_vs_fills_inv( indxs );
                dir_high                = highest_dir_per_ratio;
                
            elseif strcmp( input_wg_type, 'full' ) == true
                % normal design
                
                % first narrow down the space to take the datapoints with the
                % maximum directivity per bottom fill factor (so for each row
                % on my design space)

                indx_full_space   = obj.fill_top_bot_ratio > invert_normal_top_bot_ratio_thresh - 0.01;

                % grab variables (remember dimensions are bot fill x top bot
                % ratio)
                topbot_ratio_vs_fills_full   = topbot_ratio_mesh( :, indx_full_space );
                bot_fills_vs_fills_full      = bot_fills_mesh( :, indx_full_space );
                directivities_vs_fills_full  = obj.directivities_vs_fills( :, indx_full_space );
                angles_vs_fills_full         = obj.angles_vs_fills( :, indx_full_space );      
                periods_vs_fills_full        = obj.periods_vs_fills( :, indx_full_space );
                offsets_vs_fills_full        = obj.offsets_vs_fills( :, indx_full_space );
                scatter_str_vs_fills_full    = obj.scatter_str_vs_fills( :, indx_full_space );
                k_vs_fills_full              = obj.k_vs_fills( :, indx_full_space );
                GC_vs_fills_full             = obj.GC_vs_fills( :, indx_full_space );
                
                % remove datapoints where the angle deviates beyond some
                % angle, by artifically setting the directivity to be very
                % low
                directivities_vs_fills_full( abs( angles_vs_fills_full - obj.optimal_angle ) > 0.5 ) = 1e-6;
                
                % for each bottom fill, pick the datapoint with the highest
                % directivity
                [ highest_dir_per_bot, indx_highest_dir_per_bot ] = max( directivities_vs_fills_full, [], 2 );
                % gotta use linear indexing
                indxs                   = sub2ind( size(angles_vs_fills_full), 1:size(angles_vs_fills_full,1), indx_highest_dir_per_bot.' );
                angles_high_dir         = angles_vs_fills_full( indxs );
                periods_high_dir        = periods_vs_fills_full( indxs );
                offsets_high_dir        = offsets_vs_fills_full( indxs );
                scatter_strs_high_dir   = scatter_str_vs_fills_full( indxs );
                k_high_dir              = k_vs_fills_full( indxs );
                topbot_ratio_high_dir   = topbot_ratio_vs_fills_full( indxs );
                bot_fills_high_dir      = bot_fills_vs_fills_full( indxs );
                GC_high_dir             = GC_vs_fills_full( indxs );
                dir_high                = highest_dir_per_bot;
                
            end     % end if strcmp( input_wg_type, 'bottom' )
            
            % DEBUG plot the resulting picked out datapoints
            % angles
            figure;
            plot( 1:length(angles_high_dir), angles_high_dir, '-o' );
            title('chosen datapoints, angles');
            makeFigureNice();
            % periods
            figure;
            plot( 1:length(periods_high_dir), periods_high_dir, '-o' );
            title('chosen datapoints, periods');
            makeFigureNice();
            % offsets
            figure;
            plot( 1:length(offsets_high_dir), offsets_high_dir, '-o' );
            title('chosen datapoints, offsets');
            makeFigureNice();
            % scattering strengths
            figure;
            plot( 1:length(scatter_strs_high_dir), scatter_strs_high_dir, '-o' );
            title('chosen datapoints, scattering strengths');
            makeFigureNice();
            % k real
            figure;
            plot( 1:length(k_high_dir), real(k_high_dir), '-o' );
            title('chosen datapoints, k real');
            makeFigureNice();
            % k imag
            figure;
            plot( 1:length(k_high_dir), imag(k_high_dir), '-o' );
            title('chosen datapoints, k imaginary');
            makeFigureNice();
            % top bottom ratio
            figure;
            plot( 1:length(topbot_ratio_high_dir), topbot_ratio_high_dir, '-o' );
            title('chosen datapoints, top/bottom ratio');
            makeFigureNice();
            % bottom fill
            figure;
            plot( 1:length(bot_fills_high_dir), bot_fills_high_dir, '-o' );
            title('chosen datapoints, bottom fill');
            makeFigureNice();
            % directivity
            figure;
            plot( 1:length(dir_high), 10*log10(dir_high), '-o' );
            title('chosen datapoints, directivity (dB)');
            makeFigureNice();
            
%             % DEBUG plot bot fills vs topbot ratio
%             figure;
%             plot( topbot_ratio_high_dir, bot_fills_high_dir, '-o' );
%             xlabel('top/bottom ratio'); ylabel('bottom fill');
%             title('DEBUG bottom fill vs top/bottom ratio');
%             makeFigureNice();
%             
%             % DEBUG plot bot fills vs top fills
%             figure;
%             plot( bot_fills_high_dir, topbot_ratio_high_dir .* bot_fills_high_dir, '-o' );
%             xlabel('bottom fill'); ylabel('top fill');
%             title('DEBUG bottom fill vs top/bottom ratio');
%             makeFigureNice();
%             
%             p = polyfit( topbot_ratio_high_dir, bot_fills_high_dir, 2 )
%             
%             x = 0:0.01:1;
%             y = -(x.^3)/2 + x;
%             figure;
%             plot(x, y);
%             xlabel('bot'); ylabel('top');
            
            % now match these data points to the desired alpha
            % starting point
            start_alpha_des     = 1e-5;
            [~, indx_x]         = min(abs( alpha_des - start_alpha_des ) );
            cur_x               = xvec(indx_x);
            
            % final synthesized variables
            obj.dir_synth                   = [];
            obj.bot_fill_synth              = [];
            obj.top_bot_fill_ratio_synth    = [];
            obj.period_synth                = [];
            obj.offset_synth                = [];
            obj.angles_synth                = [];
            obj.scatter_str_synth           = [];
            obj.k_synth                     = [];
            obj.GC_synth                    = {};
            obj.des_scatter_synth           = [];
            
            % flag for switching to using max scattering strength
            saturate_scatter_str_to_max = false;
 
            ii = 1;
            while cur_x < xvec(end)
                % build grating one cell at a time
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter                 = alpha_des(indx_x);                                        % desired alpha
                if des_scatter  > max( scatter_strs_high_dir )
                    % desired scattering strength too high, gotta saturate
                    saturate_scatter_str_to_max = true;
                end
                if ~saturate_scatter_str_to_max
                    [~, indx_closest_scatter]   = min( abs(scatter_strs_high_dir - des_scatter) );          % index of closest scatter design 
                else
                    [~, indx_closest_scatter]   = max( scatter_strs_high_dir );                             % saturate to max
                end
                
                % save parameters
                obj.dir_synth(ii)                   = dir_high( indx_closest_scatter );
                obj.bot_fill_synth(ii)              = bot_fills_high_dir( indx_closest_scatter );
                obj.top_bot_fill_ratio_synth(ii)    = topbot_ratio_high_dir( indx_closest_scatter );
                obj.offset_synth(ii)                = offsets_high_dir( indx_closest_scatter );
                obj.period_synth(ii)                = periods_high_dir( indx_closest_scatter );
                obj.angles_synth(ii)                = angles_high_dir( indx_closest_scatter );
                obj.scatter_str_synth(ii)           = scatter_strs_high_dir( indx_closest_scatter );
                obj.k_synth(ii)                     = k_high_dir( indx_closest_scatter );
                obj.GC_synth{ii}                    = GC_high_dir{ indx_closest_scatter };
                obj.des_scatter_synth(ii)           = des_scatter;
                
                % move onto next
                cur_x       = cur_x + obj.period_synth(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                ii          = ii + 1;
                
            end     % end for ii = 1:ncells
            
            
        end     % end function generateFinalDesignGaussian()
        
        
        function obj = runFinalDesignEME(obj, MFD)
            % runs the final design's index distribution in EME and saves
            % some coupling parameters
            %
            % Inputs:
            %   obj
            %       takes the synthesized params as input (meaning this function
            %       can only be run after a synthesis run)
            %   MFD
            %       mode field diameter of fiber, in units 'units'
            %
            % sets these object fields:
            %   final_index
            %   obj.final_design.max_coupling_angle
            %   obj.final_design.max_coupling_offset    
            %   obj.final_design.power_reflection  
            %   obj.final_design.eme_obj 
            %   obj.final_design.final_index
            %
            % also opens the EME UI
            
            % input waveguide length
            in_wg_len   = obj.discretization * 50;                          % input waveguide length
            
            % number of cells
            n_cells = length(obj.dir_synth);
            
            
            % make input waveguide
            if strcmp( obj.input_wg_type, 'bottom' )
                % inverted design, bottom waveguide input only
                input_waveguide = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                                in_wg_len, ...
                                                0, ...
                                                1, ...
                                                0 );
            elseif strcmp( obj.input_wg_type, 'full' )
                % normal design, thick waveguide input
                input_waveguide = obj.h_makeGratingCell(  obj.convertObjToStruct(), ...
                                                in_wg_len, ...
                                                1, ...
                                                1, ...
                                                0 );
            end
            
            % lets stitch together the index distribution
            % i'm curious to see what it looks like
            obj.final_index = input_waveguide.N;
%             obj.final_index = obj.GC_synth{1}.N;
            for ii = 1:n_cells
               
                obj.final_index = [ obj.final_index, obj.GC_synth{ii}.N ];
                
            end


            % Set Up Simulation
            % note that emeSim uses 'z' as propagation direction and 'x'
            % as transverse (synthGrating uses 'x' and 'y' respectively)
            % and units are in um
            um              = 1e6;
            disc_eme        = obj.discretization * obj.units.scale * um .* [ 1, 1 ];        % [x,z]
            pol             = 0;                                                            % 0 for TE, 1 for TM
            xf              = obj.domain_size(1) * obj.units.scale * um;                    % in um (transverse domain)
            zf              = size(obj.final_index,2) * disc_eme(2);                        % in um (longitudinal domain)
            lambda_um       = obj.lambda * obj.units.scale * um;                            % wl in um
            eme_obj         = emeSim(   'discretization', disc_eme, ...
                                        'pml', 0.2, ...
                                        'domain', [xf, zf], ...
                                        'backgroundIndex', obj.background_index, ...
                                        'wavelengthSpectrum', [lambda_um lambda_um 0.1], ...
                                        'debug', 'no',...                   
                                        'polarization', pol );
                                        
            % replace the dielectric in the eme object
            eme_obj.diel = obj.final_index;
            
            % run EME sim
            % Converts the dielectric distribution into layers for eigen mode expansion
            eme_obj = eme_obj.convertDiel();   
            % Runs simulation
            eme_obj = eme_obj.runSimulation('plotSource','yes');      
            % compute fiber overlap
            z_offset    = 0 : 0.25 : zf;
            angle_vec   = obj.optimal_angle - 15 : 0.25 : obj.optimal_angle + 15;
            eme_obj     = eme_obj.fiberOverlap( 'zOffset', z_offset,...
                                                'angleVec', angle_vec,...
                                                'MFD', MFD * obj.units.scale * um,...
                                                'overlapDir', obj.coupling_direction, ...
                                                'nClad', obj.background_index );
                                        
            % DEBUG show results
            gratingUI(eme_obj);
            
            % save final results
            final_design.max_coupling_angle     = eme_obj.fiberCoup.optAngle;
            final_design.max_coupling_offset    = eme_obj.fiberCoup.optZOffset;
            final_design.max_coupling_eff       = eme_obj.fiberCoup.optCoup;
            final_design.power_reflection       = eme_obj.scatterProperties.PowerRefl(1,1);
            final_design.eme_obj                = eme_obj;
            final_design.final_index            = obj.final_index;
            final_design.input_wg_type          = obj.input_wg_type;
            final_design.MFD                    = MFD;
            obj.final_design                    = final_design;
            
            % calculate final up/down directivity
            obj = obj.calc_final_design_directivity();
            
        end     % end runFinalDesignEME()
        
        
        function obj = runFinalDesignEME_fiber_overlap(obj, MFD, fiber_offsets, angles)
            % mostly for debugging purposes
            % only runs the fiber oveerlap and not the EME simulation
            % because its faster lol
            %
            % Inputs
            %   MFD
            %       units of 'units'
            %   fiber_offsets
            %       units of 'units'
            %   angles
            %       deg.
            
            eme_obj = obj.final_design.eme_obj;
            
            % compute fiber overlap
            um      = 1e6;
            eme_obj = eme_obj.fiberOverlap( 'zOffset', fiber_offsets * obj.units.scale * um,...
                                            'angleVec', angles,...
                                            'MFD', MFD * obj.units.scale * um,...
                                            'overlapDir', obj.coupling_direction, ...
                                            'nClad', obj.background_index );
                                        
            % save final results
            final_design.max_coupling_angle     = eme_obj.fiberCoup.optAngle;
            final_design.max_coupling_offset    = eme_obj.fiberCoup.optZOffset;
            final_design.max_coupling_eff       = eme_obj.fiberCoup.optCoup;
            final_design.power_reflection       = eme_obj.scatterProperties.PowerRefl(1,1);
            final_design.eme_obj                = eme_obj;
            final_design.final_index            = obj.final_index;
            obj.final_design                    = final_design;
            
        end     % end runFinalDesignEME_fiber_overlap()
        
        
        function [obj, up_down_directivity] = calc_final_design_directivity(obj)
            % mostly for debugging purposes
            % calculates the directivity of the final designed grating,
            % using the field from the simulated eme obj
            
            % repeating code from twoLevelGratingCell which I don't like to
            % do this is bad prog. practice but whatever
            
            % define constants (all eme obj units are in um)
            mu0     = 4*pi*1e-7;                % units of H/m
            mu0     = mu0 * obj.units.scale;    % units of H/(units)
            c       = 3e8;                      % units of m/s
            c       = c/obj.units.scale;        % units of (units)/s
            omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
       
            % grab field and coordinates
            % defined as Ey in eme obj but Ez in my def.
            Ez = obj.final_design.eme_obj.fullFields.Ey;
            dy = obj.final_design.eme_obj.domain.discretization(1);
            dx = obj.final_design.eme_obj.domain.discretization(2);
            
            
%             % Calc propagating power vs. x
%             
%             % H y in , on second to end-1 steps, using entire E_z
%             % dimensions are y (transverse) vs. x
%             H_y_in = (1i/(omega0*mu0)) * ( Ez_onecell( :, 3:end ) - Ez_onecell( :, 1:end-2 ) )/(2*obj.dx);   % dx, term staggered 1
%             
%             % poynting vector in, dimensions are y vs x(2:end-1)
%             Sx = real( -1 * conj( H_y_in ) .* Ez_onecell( :, 2:end-1 ) );              % Sx = real( -Ez Hy* )
% 
%             % power per x slice
%             P_per_x_slice = sum( Sx, 1 )*obj.dy;
            

            % Calculate power radiated up
            
            % calculate H_x
            % dimensions H_x vs. y vs x
            H_x  = 1/(1i*omega0*mu0) .* ( Ez( 3:end,:) - Ez( 1:end-2,:) )/(2*dy);

            % power flowing across y
            Sy                  = real( conj( H_x ) .* Ez( 2:end-1,: ) );           % Sy = real( Ez Hx* )
            P_per_y_slice       = -sum( Sy, 2 )*dx;                                 % using Sy compmonent, for some reason up/down are flipped

            [~, indx_max_up]    = max( P_per_y_slice );
            [~, indx_max_down]  = min( P_per_y_slice );
            Sy_up               = Sy( indx_max_up, : );                                 % dimensions Sy vs. x
            Sy_down             = Sy( indx_max_down, : );
            
%             % save to obj
% %             obj.Sx              = Sx;
% %             obj.P_per_x_slice 	= P_per_x_slice;
%             obj.Sy              = Sy;
%             obj.P_per_y_slice   = P_per_y_slice;

            % calc radiated power (Sy only, not including angle of
            % incidence)
            P_rad_up    = sum( abs( Sy_up(:) ) ) * dx;
            P_rad_down  = sum( abs( Sy_down(:) ) ) * dx;
            
            % calc directivity save to object
            up_down_directivity                  = P_rad_up/P_rad_down;
            obj.final_design.up_down_directivity = up_down_directivity;
            
        end     % end calc_final_design_directivity()
        
        
        function [obj, max_eff_vs_MFD] = calc_eff_vs_MFD(obj, MFDs)
            % mostly for debugging
            % plotting max efficiency vs. different MFDs
            %
            % returns max efficiency vs. mfd
            %
            % Inputs
            %   MFDs
            %       type: double, vector
            %       desc: mode field diameters to try, units 'units'
            %
            % Outputs
            %   max_eff_vs_MFD
            %       type: 3d tensor, double
            %       desc: eff vs. 
            
            % init vars
            max_eff_vs_MFD  = zeros(size(MFDs));
            fiber_offsets   = ( 0:0.5:12 ) * 1e-6 / obj.units.scale;        % units 'units'
            angles          = 0:0.1:45;

            for ii = 1:length(MFDs)

                fprintf('mfd iteration %i of %i\n', ii, length(MFDs) );

                % run EME sim
                obj = obj.runFinalDesignEME_fiber_overlap( MFDs(ii), fiber_offsets, angles );
                % save result
                max_eff_vs_MFD(ii) = obj.final_design.max_coupling_eff;

            end
            
        end     % end calc_eff_vs_MFD()
        
        
        function [obj, eff_vs_wavelength_angle] = calc_eff_vs_wavelength(obj, wavelengths)
            % mostly for debugging
            % plotting eff. vs. wavelength and angle
            %
            % i dont really want to do this right now haha
            
            % init vars
            angles                  = 0:0.1:45;
            eff_vs_wavelength_angle = zeros( length(wavelengths), length(angles) );     % dimensions eff. vs. wl vs. angle
            
            for ii = 1:length(wavelengths)
               
                fprintf('wavelength iteration %i of %i\n', ii, length(wavelengths) );
                
                % run EME sim
                obj = obj.runFinalDesignEME_fiber_overlap( MFDs(ii), fiber_offsets, angles );
                % save result
                max_eff_vs_MFD(ii) = obj.final_design.max_coupling_eff;
                
            end
            
            
        end     % end calc_eff_vs_wavelength()
        
        
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
% currently deprecated
% makes and returns a c_twoLevelGratingCell object
% 
% inputs:
%   synth_obj
%       type: c_synthGrating object AS STRUCT
%       desc: c_synthGrating object AS STRUCT
%   period
%       type: double, scalar
%       desc: period of the grating cell
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

                        
% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );
                        
% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = 100.0 / ( 1e9 * synth_obj.units.scale ) * ones(1,2);      % hardcoded to 100nm
wg_index        = [3.45, 3.45];                                             % hardcoded to 3.45
wg_min_y        = [ domain_size(1)/2, domain_size(1)/2-wg_thick(1) ];
wgs_duty_cycles = [ fill_top, fill_bot ];
wgs_offsets     = [ 0, offset_ratio*period ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, wg_index, ...
                                        wgs_duty_cycles, wgs_offsets );
            
end




















































