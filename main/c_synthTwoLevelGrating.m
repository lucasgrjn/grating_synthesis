classdef c_synthTwoLevelGrating < c_synthGrating
% Synthesizes a high directivity 2-level grating in an arbitrary process
%
% Authors: bohan zhang
%
% this documentation needs to be updated
%
% Prerequisites/dependencies
%   - c_synthGrating.m
%   - c_twoLevelGratingCell.m
%   - the utility folder
%
%   The user should define their own custom grating unit cell
%   drawing function.
%   HOWEVER, this function MUST have the following inputs and outputs, IN
%   ORDER:
%       function GC = your_makeGratingCell_function( dxy, units, lambda, background_index, y_domain_size, ...
%                                          period, fill_top, fill_bot,
%                                          offset_ratio );
%           % makes and returns a c_twoLevelGratingCell object
%           inputs(outdated):
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
%   'data_notes'
%       type: string
%       desc: optional verbose notes/descriptor for this simulation
%
%   'h_makeGratingCell'
%       type: function handle
%       desc: handle to grating drawing function
%
% Examples:
%
%   % make synthesis object
%   synth_obj = c_synthTwoLevelGrating(   'discretization',    disc, ...
%                                       'units',             units,   ...
%                                       'lambda',            lambda, ...
%                                       'background_index',  background_index,    ...
%                                       'y_domain_size',     y_domain_size, ...
%                                       'optimal_angle',     optimal_angle, ...
%                                       'data_notes',        data_notes, ...
%                                       'coupling_direction', coupling_direction, ...
%                                       'h_makeGratingCell', @f_makeGratingCell_45RFSOI ...
%                                       );



    properties
        
        % inherits the following from c_synthGrating
        % discretization;         % dx and dy
        % units;                  % units, verbose, 'm' or 'mm', or 'um', or 'nm'
        %                         % has fields 'name' and 'scale'
        % lambda;                 % center wavelength
        % background_index;       % background index
        % y_domain_size;          % transverse domain size (y coordinate aka vertical dimension)
        % inputs;                 % saves input settings for user reference
        % start_time;             % time when object was created, 'YEAR-month-day hour-min-sec'
        % data_notes;             % verbose notes of what current sweep is doing     
        % h_makeGratingCell;      % handle to the grating cell making function
        % coupling_direction;     % either 'up' or 'down
        % optimal_angle;          % angle to optimize for, deviation from the normal, in deg.
        % sweep_variables;
        % synthesized_design;

        % fields saved in sweep_variables:
        % dimensions are top/bot ratio vs. bot fill
        % fill_tops
        % fill_bots
        % fill_top_bot_ratio
        % directivities_vs_fills
        % angles_vs_fills
        % scatter_str_vs_fills
        % periods_vs_fills
        % offsets_vs_fills - are these... ratios? or absolute values?
        % k_vs_fills
        % GC_vs_fills
        % dir_b4_period_vs_fills
                            
        % fields saved in synthesized_design:
        % dimensions are vs. cell #
        % dir                 
        % bot_fill             
        % top_bot_fill_ratio  
        % top_fill             
        % period               
        % offset        % absolute offset    
        % angles              
        % scatter_str         
        % k                 
        % GC               
        % des_scatter  
        % N
        % x_coords
        % y_coords
        % input_wg_type
        % use_min_feat_size     currently either true or false, doesn't save values
        
        % debug options
        % a struct that holds debugging options
        % currently saved fields:
        %   verbose
        debug_options;
        
%         input_wg_type;  % 'bottom' or 'full

        chosen_cells;
        
    end
    
    methods
        
        function obj = c_synthTwoLevelGrating(varargin)
            % Constructor
            % See top comments for input documentation
            
            % call c_synthGrating constructor
            obj = obj@c_synthGrating(varargin{:});
            
        end     % end constructor()

        

        % -----------------
        % Function generate_design_space()
        %
        %> Generates the design space vs. fill factors
        %> picking optimum offset and period for highest directivity and
        %> closest angle to desired
        %>
        %> This version runs versus top and bottom fill factors
        %>
        %> based on synthesizeGaussianGrating
        %>
        %> would be good to implement: option to save GC data or not
        %>
        %> Inputs:
        %>   fill_bots
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing
        %>   fill_tops
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing
        %>   verbose
        %>       type: double, array
        %>       desc: OPTIONAL Currently mostly for testing, spits out
        %>             a bunch of stuff to the prompt
        function obj = generate_design_space( obj, fill_bots, fill_tops, verbose )
            
            tic;
            fprintf('Sweeping fill factors for directivity and angle...\n');
            
            % set verbose options
            if ~exist('verbose', 'var')
                % default verbose off
                obj.debug_options.verbose = false;
            else
                obj.debug_options.verbose = verbose;
            end
            
            % set fill factors and offsets
            if ~exist('fill_bots', 'var')
                % default fill
                fill_bots   = fliplr( 0.025:0.025:0.975 );
            end
            if ~exist('fill_tops', 'var')
                % default fill
                fill_tops   = fliplr( 0.025:0.025:0.975 );
            end
            fill_top_bot_ratio  = [];
            guess_offset        = 0;
            
            % number of times to iterate over optimizatio loop
            n_optimize_loops = 2;
            
            % sort fills so they are in descending order
            fill_bots = sort( fill_bots, 'descend' );
            fill_tops = sort( fill_tops, 'descend' );
            
            % save fills and offsets
            obj.sweep_variables.fill_tops           = fill_tops;
            obj.sweep_variables.fill_bots           = fill_bots;
            obj.sweep_variables.fill_top_bot_ratio  = fill_top_bot_ratio;
            
            % initialize saving variables
            directivities_vs_fills  = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            angles_vs_fills         = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            periods_vs_fills        = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            offsets_vs_fills        = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill, this is offset ratio
            scatter_str_vs_fills    = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            k_vs_fills              = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            GC_vs_fills             = cell( length( fill_bots ), length( fill_tops ) );      % dimensions bot fill vs. top fill
            dir_b4_period_vs_fills  = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            prad_pin_vs_fills       = zeros( length( fill_bots ), length( fill_tops ) );     % dimensions bot fill vs. top fill
            
            % make grating cell, assuming both layers are filled
            waveguide = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               2*obj.discretization, ...
                                               1.0, 1.0, 0.0 );
            
            % run waveguide simulation
            % sim settings
            guess_n             = 1.0 * max( waveguide.N(:) );                                      % guess index. I wonder if there's a better guessk for this?
            guessk              = guess_n * 2*pi/obj.lambda;                                        % units rad/'units'
            num_wg_modes        = 5;
            BC                  = 0;                                                                % 0 = PEC
            pml_options_wg      = [0, 200, 20, 2];                                                  % now that I think about it... there's no reason for the user to set the pml options
            % run sim
            waveguide   = waveguide.runSimulation( num_wg_modes, BC, pml_options_wg, obj.k0, guessk );
            
            % update guessk (units rad/'units')
            guessk = waveguide.k;
            
            % grab waveguide k
            waveguide_k = waveguide.k;                                                              % units of rad/'units'
            
            % calculate analytical period which would approximately phase
            % match to desired output angle
            k0              = obj.background_index * ( 2*pi/obj.lambda );
            kx              = k0 * sin( (pi/180) * obj.optimal_angle );
            guess_period    = 2*pi/(waveguide_k- kx);                                               % units of 'units'
            
            % snap period to discretization
            guess_period = obj.discretization * round(guess_period/obj.discretization);
            
            % ugh this is really annoying but i have to - extend the
            % waveguide's e z overlap
            [ waveguide, e_z_overlap_ext ]  = ...
                waveguide.stitch_E_field( waveguide.Phi, real(waveguide.k), round(guess_period/waveguide.domain_size(2)) );
            waveguide.E_z_for_overlap       = e_z_overlap_ext;
            
            % set grating solver settings
            num_modes   = 5;
            BC          = 0;                                                % 0 = PEC
            pml_options = [1, 100, 20, 2]; 
            OPTS        = struct( 'mode_to_overlap', e_z_overlap_ext );
            sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);
            
            % initially start with waveguide GC
            guess_GC = waveguide;
                
            % first fill in the right side of the domain
            for i_ff_bot = 1:length( fill_bots )

                % print iteration
                fprintf('Right side of domain, fill iteration %i of %i\n', i_ff_bot, length(fill_bots) );

                % Optimize period and offset
                if fill_bots(i_ff_bot) < 1 && fill_tops(1) < 1
                    % Only run optimization if theres a perturbation of
                    % both layers
                    for ii = 1:n_optimize_loops
                        % run optimization loop
                        [ ~, ...
                          guess_period, ...
                          guess_offset, ...
                          directivities_vs_fills( i_ff_bot, 1 ), ...
                          angles_vs_fills( i_ff_bot, 1 ), ...
                          scatter_str_vs_fills( i_ff_bot, 1 ), ...
                          guess_GC, ...
                          guessk, ...
                          dir_b4_period_vs_fills( i_ff_bot, 1 ) ...
                          ] = ...
                            obj.optimize_period_offset( guess_offset, ...
                                                      fill_tops(1), ...
                                                      fill_bots(i_ff_bot), ...
                                                      guess_period,...
                                                      guessk, ...
                                                      sim_opts, ...
                                                      guess_GC );
                    end

                    % save the optimized period, offset, GC, and k
                    periods_vs_fills( i_ff_bot, 1 ) = guess_period;
                    offsets_vs_fills( i_ff_bot, 1 ) = guess_offset;
                    GC_vs_fills{ i_ff_bot, 1 }      = guess_GC;
                    k_vs_fills( i_ff_bot, 1 )       = guessk;
                    
                    if strcmp( obj.coupling_direction, 'up' )
                        prad_pin_vs_fills( i_ff_bot, 1 ) = guess_GC.P_rad_up/guess_GC.P_in;
                    else
                        prad_pin_vs_fills( i_ff_bot, 1 ) = guess_GC.P_rad_down/guess_GC.P_in;
                    end
                    
                else
                    % at least one of the layers is not perturbed, so there is no optimization to run
                    % save dummy values
                    directivities_vs_fills( i_ff_bot, 1 )   = 1;
                    angles_vs_fills( i_ff_bot, 1 )          = 0;
                    scatter_str_vs_fills( i_ff_bot, 1 )     = 0;
                    periods_vs_fills( i_ff_bot, 1 )         = guess_period;
                    offsets_vs_fills( i_ff_bot, 1 )         = 0;
                    k_vs_fills( i_ff_bot, 1 )               = guessk;
                    GC_vs_fills{ i_ff_bot, 1 }              = waveguide;
                    dir_b4_period_vs_fills( i_ff_bot, 1 )   = 1;
                    prad_pin_vs_fills( i_ff_bot, 1 )        = 0;

                end     % end if fill_top < 1

                toc;

            end     % end initial domain sweep
            
            
            % for parallel processing, grab these variables before entering
            % parfor
            offsets_vs_fills_1 = offsets_vs_fills(:,1);
            periods_vs_fills_1 = periods_vs_fills(:,1);
            k_vs_fills_1       = k_vs_fills(:,1);
            GC_vs_fills_1      = GC_vs_fills(:,1);
            % calc number of loops, also necessary apparently for parfor
            n_fill_bots    = length(fill_bots);
            n_fill_tops    = length(fill_tops);

            % start up parallel pool
            obj = obj.start_parpool();
            
            % now fill in the rest of the domain
            parfor i_ff_bot = 1:n_fill_bots
                % For each bottom fill factor
    
                fprintf('Main parfor iteration %i of %i\n', i_ff_bot, n_fill_bots);
                
                % grab starting guess period and k
                guess_period    = periods_vs_fills_1( i_ff_bot );
                guessk          = k_vs_fills_1( i_ff_bot );
                guess_GC        = GC_vs_fills_1{ i_ff_bot };
                guess_offset    = offsets_vs_fills_1( i_ff_bot );

                % grab bottom fill
                fill_bot = fill_bots( i_ff_bot );

                for i_ff_top = 2:n_fill_tops
                    % for each top/bottom fill factor ratio

                    fprintf('Fill factor ratio %i of %i, main parfor iteration %i of %i\n', i_ff_top, n_fill_tops, i_ff_bot, n_fill_bots);
                        
                    fill_top = fill_tops(i_ff_top);
                    
                    % Optimize period and offset
%                     fill_top = fill_top_bot_ratio(i_ff_ratio) * fill_bot;
                    if fill_top < 1
                        % Only run optimization if theres a perturbation
                        
                        for ii = 1:n_optimize_loops
                            % run optimization loop
                            [ ~, ...
                              guess_period, ...
                              guess_offset, ...
                              best_directivity, ...
                              best_angle, ...
                              best_scatter_str, ...
                              guess_GC, ...
                              guessk, ...
                              dir_b4_period_vs_fill ...
                              ] = ...
                                obj.optimize_period_offset( guess_offset, ...
                                                          fill_top, ...
                                                          fill_bot, ...
                                                          guess_period,...
                                                          guessk, ...
                                                          sim_opts, ...
                                                          guess_GC );
                        end
                            
                        % save data
                        periods_vs_fills( i_ff_bot, i_ff_top )        = guess_period;
                        offsets_vs_fills( i_ff_bot, i_ff_top )        = guess_offset;
                        directivities_vs_fills( i_ff_bot, i_ff_top )  = best_directivity;
                        angles_vs_fills( i_ff_bot, i_ff_top )         = best_angle;
                        scatter_str_vs_fills( i_ff_bot, i_ff_top )    = best_scatter_str;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }             = best_GC;
                        k_vs_fills( i_ff_bot, i_ff_top )              = guessk;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_top )  = dir_b4_period_vs_fill;
                        if strcmp( obj.coupling_direction, 'up' )
                            prad_pin_vs_fills( i_ff_bot, i_ff_top ) = guess_GC.P_rad_up/guess_GC.P_in;
                        else
                            prad_pin_vs_fills( i_ff_bot, i_ff_top ) = guess_GC.P_rad_down/guess_GC.P_in;
                        end
                            
                    else
                        % at least one of the layers is not perturbed, so there is no optimization to run
                        % save dummy values
                        directivities_vs_fills( i_ff_bot, i_ff_top )    = 1;
                        angles_vs_fills( i_ff_bot, i_ff_top )           = 0;
                        scatter_str_vs_fills( i_ff_bot, i_ff_top )      = 0;
                        periods_vs_fills( i_ff_bot, i_ff_top )          = guess_period;
                        offsets_vs_fills( i_ff_bot, i_ff_top )          = 0;
                        k_vs_fills( i_ff_bot, i_ff_top )                = guessk;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }               = waveguide;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_top )    = 1;
                        prad_pin_vs_fills( i_ff_bot, i_ff_top )         = 0;

                    end     % end if fill_top < 1
                end     % end for i_ff_top = ...         
            end     % end parfor i_ff_bot = ...
            
            % save variables to object
            obj.sweep_variables.directivities_vs_fills  = directivities_vs_fills;
            obj.sweep_variables.angles_vs_fills         = angles_vs_fills;
            obj.sweep_variables.scatter_str_vs_fills    = scatter_str_vs_fills;
            obj.sweep_variables.periods_vs_fills        = periods_vs_fills;
            obj.sweep_variables.offsets_vs_fills        = offsets_vs_fills;
            obj.sweep_variables.k_vs_fills              = k_vs_fills;
            obj.sweep_variables.dir_b4_period_vs_fills  = dir_b4_period_vs_fills;
            obj.sweep_variables.prad_pin_vs_fills       = prad_pin_vs_fills;

            fprintf('Done generating design space\n');
            toc;
            
        end     % end generate_design_space()
        
        
        
        function obj = start_parpool( obj )
            % convenience function for starting up a parallel pool
            
            % start a parallel pool session
            my_cluster = parcluster('local');                           % cores on compute node are "local"
            if getenv('ENVIRONMENT')                                    % true if this is a batch job
                my_cluster.JobStorageLocation = getenv('TMPDIR');       % points to TMPDIR
            end

            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if ~isempty(poolobj)
                % shut down previously made parallel pool
                delete(gcp('nocreate'));
            end
            % Automatically selects number of works to the max that the
            % cluster supports
            parpool(my_cluster, my_cluster.NumWorkers);
            
        end     % end start_parpool()

        
        
        function [  obj, best_period, best_offset, best_directivity, best_angle, best_scatter_str, ...
                    best_GC, best_k, dir_b4_period_vs_fill, DEBUG ] ...
                    = optimize_period_offset(obj, guess_offset, fill_top, fill_bot, guess_period, guessk, sim_opts, guess_gc )
            % for given fill ratios, optimizes period and offset for best angle/directivity
            %
            % Inputs:
            %   guess_offset
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
            %       type: scalar, double
            %       desc: optimal period
            %   best_offset
            %       type: scalar, double
            %       desc: absolute value of offset, in 'units'
            %   best_directivity
            %       type: scalar, double
            %       desc: optimal directivity
            %   best_angle
            %       type: scalar, double
            %       desc: optimal angle
            %   best_scatter_str
            %       type: scalar, double
            %       desc: optimal scatering strength
            %   best_GC
            %       type: c_twoLevelGratingCell
            %       desc: two level grating cell
            %   best_k
            %       type: scalar, double
            %       desc: optimal k
            %   dir_b4_period_vs_fill
            %       type: scalar, double
            %       desc: directivity b4 period sweep, mostly for debugging purposes 
            %   DEBUG
            %       type: struct
            %       desc: as name suggests, holds fields for debugging
            
            directivities = [];
            k_vs_offset   = [];
            angles        = [];
            GC_vs_offset  = {};

            % grab mode to overlap with
            OPTS = struct( 'mode_to_overlap', guess_gc.E_z_for_overlap );

            % new version of sweepign offsets, local search
            cur_offset   = guess_offset;
            offsets      = [];
            delta_offset = obj.discretization;  % start with positive delta offset
            i_offset     = 1;
            while true
               
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               guess_period, ...
                                               fill_top, ...
                                               fill_bot, ...
                                               cur_offset./guess_period );
                
                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, obj.k0, guessk, OPTS );
                
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
                
                offsets( i_offset ) = cur_offset;

                % Check if we need to switch directions
                if i_offset > 1
                    
                    if delta_offset > 0 && directivities( i_offset ) < directivities( i_offset-1 )
                        % Switch directions, we hit a max.
                        delta_offset            = -delta_offset;
                        OPTS.mode_to_overlap    = guess_gc.E_z_for_overlap;
                        cur_offset              = guess_offset;
                    elseif delta_offset < 0 && directivities( i_offset ) < max( directivities( 1:end-1) )
                        % Quit loop, we're done
                        break;
                    end
                    
                end
                
                % update for next loop
                cur_offset      = cur_offset + delta_offset;
                i_offset        = i_offset + 1;
                    
            end

            % DEBUG field
            DEBUG.GC_vs_offset = GC_vs_offset;

            % pick best offset, to feed into period sweep loop
            [ ~, indx_best_offset ]     = max( directivities );
%             best_offset_ratio           = offset_ratios( indx_best_offset );
            best_offset_ratio           = offsets( indx_best_offset )./guess_period;      % newer local search ver.
            best_offset_k               = k_vs_offset( indx_best_offset );
            best_offset_angle           = angles( indx_best_offset );
            
            % update mode overlap
            OPTS.mode_to_overlap    = GC_vs_offset{indx_best_offset}.E_z_for_overlap;

            % here's an output variable
            dir_b4_period_vs_fill = max( directivities );

            % now sweep periods
            % decide whether to sweep larger or smaller periods
            % based on the angle
            if best_offset_angle > obj.optimal_angle
                % only sweep smaller periods
                delta_period = -obj.discretization;
            else
                % only sweep larger periods
                delta_period = obj.discretization;
            end

            % init saving variables
            angles_vs_period    = []; 
            k_vs_period         = [];
            GC_vs_period        = {};
            periods             = [];
            
            % initial period sweep values
            guessk              = best_offset_k;
            periods(1)          = guess_period;
            GC_vs_period{1}     = GC_vs_offset{indx_best_offset};
            k_vs_period(1)      = best_offset_k;
            angles_vs_period(1) = best_offset_angle;
                       
            i_period    = 2;
            while true
             
                % verbose printing
                if obj.debug_options.verbose == true
                    fprintf('Sweeping period %i\n', i_period );
                end
                
                % update period
                guess_period = guess_period + delta_period;
                
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               guess_period, ...
                                               fill_top, ...
                                               fill_bot, ...
                                               best_offset_ratio );

                % run sim
                GC = GC.runSimulation( sim_opts.num_modes, sim_opts.BC, sim_opts.pml_options, obj.k0, guessk, OPTS );

                % save angle
                if strcmp( obj.coupling_direction, 'up' )
                    % coupling direction is upwards
                    angles_vs_period( i_period ) = GC.max_angle_up;
                else
                    % coupling direction is downwards
                    angles_vs_period( i_period ) = GC.max_angle_down;
                end

                % update for next iteration
                periods(i_period)       = guess_period;
                GC_vs_period{i_period}  = GC;
                k_vs_period(i_period)   = GC.k;
                guessk                  = GC.k;
                OPTS.mode_to_overlap    = GC.E_z_for_overlap;
                
                % check for exit condition (if error in angle gets worse)
                cur_angle_err   = abs( angles_vs_period( i_period ) - obj.optimal_angle );
                prev_angle_err  = abs( angles_vs_period( i_period-1 ) - obj.optimal_angle );
                if cur_angle_err > prev_angle_err
                    % optimization over, break
                    break;
                end
                
                i_period = i_period + 1;

            end     % end period sweep

            % DEBUG field
            DEBUG.GC_vs_period = GC_vs_period;
            
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
%             best_offset = best_offset_ratio * best_period;
            best_offset = round( best_offset_ratio * best_period / obj.discretization ) .* obj.discretization;  % snap to grid
                    
        end     % end function optimize_period_offset()
              
        
        
        function obj = generate_final_design_apodized( obj, desired_field, input_wg_type, enforce_min_feat_size_func, ...
                                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % function for generating the final synthesized design
            % parameters
            %
            % uses the top and bottom fills
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter, in units 'units'
            %   input_wg_type
            %       type: string
            %       desc: 'bottom' for body only or 'full' for both layers
            %             or 'none' for no layers
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL INPUT
            %             A function that user makes which enforces min.
            %             feat size
            %             See the example function at the bottom of this
            %             file
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            %   
            %
            % Sets these fields: (not updated)
            %     obj.synthesized_design.dir                  
            %     obj.synthesized_design.bot_fill            
            %     obj.synthesized_design.top_bot_fill_ratio  
            %     obj.synthesized_design.top_fill            
            %     obj.synthesized_design.period              
            %     obj.synthesized_design.offset              
            %     obj.synthesized_design.angles             
            %     obj.synthesized_design.scatter_str        
            %     obj.synthesized_design.k                  
            %     obj.synthesized_design.GC                 
            %     obj.synthesized_design.des_scatter             
            
            % calcualte desired scattering
            [ obj, xvec, alpha_des, desired_field ] = obj.calculate_desired_scattering( desired_field );
            
            % select cells to use
            obj = obj.choose_unit_cells( input_wg_type, enforce_min_feat_size_func, ...
                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both );

                     
            % optimize start alpha
            [ obj, best_alpha_power ] = obj.optimize_start_alpha( xvec, alpha_des, obj.chosen_cells.chosen_dirs, ...
                                              obj.chosen_cells.chosen_bot_fills, obj.chosen_cells.chosen_top_fills, ...
                                              obj.chosen_cells.chosen_offsets, obj.chosen_cells.chosen_periods, ...
                                              obj.chosen_cells.chosen_angles, obj.chosen_cells.chosen_scatter_str, ...
                                              obj.chosen_cells.chosen_ks, desired_field );
                                          
            % DEBUG PLOTTING PREDICTED FIELD SHAPE
            [ ~, field_shape_prediction ] = obj.predict_overlap_for_optimization( ...
                                              xvec, alpha_des, obj.chosen_cells.chosen_dirs, ...
                                              obj.chosen_cells.chosen_bot_fills, obj.chosen_cells.chosen_top_fills, ...
                                              obj.chosen_cells.chosen_offsets, obj.chosen_cells.chosen_periods, ...
                                              obj.chosen_cells.chosen_angles, obj.chosen_cells.chosen_scatter_str, ...
                                              obj.chosen_cells.chosen_ks, desired_field, best_alpha_power );
%             figure('Name','predicted_field');
%             plot( xvec, alpha_des ); hold on;
%             plot( xvec, abs(desired_field).*max(alpha_des)./max(desired_field) );
%             plot( xvec, field_shape_prediction.*max(alpha_des)./max(field_shape_prediction) );
%             xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
%             legend('desired scattering strength', 'fiber mode', 'predicted field shape');
%             title('DEBUG predicted field shape');
%             makeFigureNice();  
                                          
            % synthesize final design
            [ obj, synthesized_des ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                 obj.chosen_cells.chosen_dirs, ...
                                 obj.chosen_cells.chosen_bot_fills, ...
                                 obj.chosen_cells.chosen_top_fills, ...
                                 obj.chosen_cells.chosen_offsets, ...
                                 obj.chosen_cells.chosen_periods, ...
                                 obj.chosen_cells.chosen_angles, ...
                                 obj.chosen_cells.chosen_scatter_str, ...
                                 obj.chosen_cells.chosen_ks, ...
                                 10.^best_alpha_power );
            obj.synthesized_design = catstruct( obj.synthesized_design, synthesized_des );  % in aux function
            % save input waveguide type
            obj.synthesized_design.input_wg_type = input_wg_type; 
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end function generate_final_design_apodized()
        
        function obj = generate_final_design_apodized_v2_20200612( obj, desired_field, input_wg_type, enforce_min_feat_size_func, ...
                                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % function for generating the final synthesized design
            % parameters
            %
            % new version created on 6/12/2020 that uses a different method
            % of matching total radiated power to apodize
            %
            % i don't recall if this really worked out in the end
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter, in units 'units'
            %   input_wg_type
            %       type: string
            %       desc: 'bottom' for body only or 'full' for both layers
            %             or 'none' for no layers
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL INPUT
            %             A function that user makes which enforces min.
            %             feat size
            %             See the example function at the bottom of this
            %             file
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            %   
            %
            % Sets these fields: (not updated)
            %     obj.synthesized_design.dir                  
            %     obj.synthesized_design.bot_fill            
            %     obj.synthesized_design.top_bot_fill_ratio  
            %     obj.synthesized_design.top_fill            
            %     obj.synthesized_design.period              
            %     obj.synthesized_design.offset              
            %     obj.synthesized_design.angles             
            %     obj.synthesized_design.scatter_str        
            %     obj.synthesized_design.k                  
            %     obj.synthesized_design.GC                 
            %     obj.synthesized_design.des_scatter             

            % select subset of domain to use
             if nargin < 5
                 [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks, rad_power_ratio ] ...
                           = obj.pick_design_quadrant( input_wg_type );
             else
                 % override the fills as user desired
                 [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks, rad_power_ratio ] ...
                           = obj.pick_design_quadrant( input_wg_type, fill_top_override, fill_bot_override );
             end
              
             
            % default to using top fill
            if nargin < 7
                fill_vs_top_bot_both = 'top';
            end
            
            % for each top fill, pick bottom fill with highest
            % directivity
            [ chosen_dirs_vs_top, indx_highest_dir_per_topfill ] = max( directivities, [], 1 );
            % gotta use linear indexing
            indxs_vs_top               = sub2ind( size(directivities), indx_highest_dir_per_topfill, 1:size(directivities,2)  );


            % for each bottom fill, pick top fill
            % with highest directivity
            % directivity
            [ chosen_dirs_vs_bot, indx_highest_dir_per_botfill ] = max( directivities, [], 2 );
            % gotta use linear indexing
            indxs_vs_bot                   = sub2ind( size(directivities), 1:size(directivities,1), indx_highest_dir_per_botfill.' );

            switch fill_vs_top_bot_both
                
                case 'top'
            
                    % use top fill
                    chosen_angles       = angles( indxs_vs_top );
                    chosen_periods      = periods( indxs_vs_top );
                    chosen_offsets      = offsets( indxs_vs_top );
                    chosen_scatter_str  = scatter_strs( indxs_vs_top );
                    chosen_ks           = ks( indxs_vs_top );
                    chosen_top_fills    = top_fills( indxs_vs_top );
                    chosen_bot_fills    = bot_fills( indxs_vs_top );
                    chosen_dirs         = chosen_dirs_vs_top;
                    chosen_rad_power_ratio = rad_power_ratio( indxs_vs_top );
                    
                case 'bottom'
                     
                    % use bot fill
                    chosen_angles       = angles( indxs_vs_bot );
                    chosen_periods      = periods( indxs_vs_bot );
                    chosen_offsets      = offsets( indxs_vs_bot);
                    chosen_scatter_str  = scatter_strs( indxs_vs_bot );
                    chosen_ks           = ks( indxs_vs_bot );
                    chosen_top_fills    = top_fills( indxs_vs_bot );
                    chosen_bot_fills    = bot_fills( indxs_vs_bot );
                    chosen_dirs         = chosen_dirs_vs_bot;
                    chosen_rad_power_ratio = rad_power_ratio( indxs_vs_bot );
                    
                case 'both'
                    
                    % append bot fill
                    chosen_angles       = [ angles( indxs_vs_top ), angles( indxs_vs_bot ) ];
                    chosen_periods      = [ periods( indxs_vs_top ), periods( indxs_vs_bot ) ];
                    chosen_offsets      = [ offsets( indxs_vs_top ), offsets( indxs_vs_bot) ];
                    chosen_scatter_str  = [ scatter_strs( indxs_vs_top ), scatter_strs( indxs_vs_bot ) ];
                    chosen_ks           = [ ks( indxs_vs_top ), ks( indxs_vs_bot ) ];
                    chosen_top_fills    = [ top_fills( indxs_vs_top ), top_fills( indxs_vs_bot ) ];
                    chosen_bot_fills    = [ bot_fills( indxs_vs_top ), bot_fills( indxs_vs_bot ) ];
                    chosen_dirs = [ chosen_dirs_vs_top, chosen_dirs_vs_bot.' ];
                    chosen_rad_power_ratio = [ rad_power_ratio( indxs_vs_top ), rad_power_ratio( indxs_vs_bot ) ]; 
            
            end
            
%             chosen_dirs = chosen_dirs_vs_top;
%             chosen_dirs = [ chosen_dirs_vs_top, chosen_dirs_vs_bot.' ];
            
            % sort on radiated power ratio
            [ chosen_rad_power_ratio, indx_sort ] = sort( chosen_rad_power_ratio );
            chosen_angles       = chosen_angles( indx_sort );
            chosen_periods      = chosen_periods( indx_sort );
            chosen_offsets      = chosen_offsets( indx_sort );
            chosen_ks           = chosen_ks( indx_sort );
            chosen_top_fills    = chosen_top_fills( indx_sort );
            chosen_bot_fills    = chosen_bot_fills( indx_sort );
            chosen_dirs         = chosen_dirs( indx_sort );
            chosen_scatter_str = chosen_scatter_str( indx_sort );  

            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;
                
                indices_to_keep = [];
                for ii = 1:length( chosen_periods )
                    if enforce_min_feat_size_func( chosen_periods(ii), chosen_top_fills(ii), chosen_bot_fills(ii), chosen_offsets(ii) ) == true
                        indices_to_keep(end+1) = ii;
                    end
                end
                
                chosen_angles       = chosen_angles(indices_to_keep);
                chosen_periods      = chosen_periods(indices_to_keep);
                chosen_offsets      = chosen_offsets(indices_to_keep);
                chosen_scatter_str  = chosen_scatter_str(indices_to_keep);
                chosen_ks           = chosen_ks(indices_to_keep);
                chosen_top_fills    = chosen_top_fills(indices_to_keep);
                chosen_bot_fills    = chosen_bot_fills(indices_to_keep);
                chosen_dirs         = chosen_dirs(indices_to_keep);
                chosen_rad_power_ratio = chosen_rad_power_ratio(indices_to_keep);
                        
            end
            
            chosen_params = struct( 'chosen_angles', chosen_angles, 'chosen_periods', chosen_periods, ...
                                'chosen_offsets', chosen_offsets, 'chosen_scatter_str', chosen_scatter_str, ...
                                'chosen_ks', chosen_ks, 'chosen_top_fills', chosen_top_fills, ...
                                'chosen_bot_fills', chosen_bot_fills, 'chosen_dirs', chosen_dirs, ...
                                'chosen_rad_power_ratio', chosen_rad_power_ratio );
                
            OPTS = struct('fig_size', [1200, 800]);
            figure('Name', 'chosen_datapoints');
            % chosen angles
            subplot(3,3,1);
            plot( 1:length(chosen_angles), chosen_angles, '-o' );
            title('chosen datapoints, angles');
            xlim( [1, length(chosen_angles) ] );
            makeFigureNice( OPTS );
            % chosen periods
            subplot(3,3,2);
            plot( 1:length(chosen_periods), chosen_periods, '-o' );
            title('chosen datapoints, periods');
            makeFigureNice( OPTS );
            % chosen offsets
            subplot(3,3,3);
            plot( 1:length(chosen_offsets), chosen_offsets, '-o' );
            title('chosen datapoints, offsets');
            makeFigureNice( OPTS );
            % chosen scatter strength
            subplot(3,3,4);
            plot( 1:length(chosen_rad_power_ratio), chosen_rad_power_ratio, '-o' );
            title('chosen datapoints, radiated power ratio');
            makeFigureNice( OPTS );
            % chosen k real
            subplot(3,3,5);
            plot( 1:length(chosen_ks), real(chosen_ks), '-o' );
            title('chosen datapoints, k real');
            makeFigureNice( OPTS );
            % chosen k imag
            subplot(3,3,6);
            plot( 1:length(chosen_ks), imag(chosen_ks), '-o' );
            title('chosen datapoints, k imaginary');
            makeFigureNice( OPTS );
            % chosen top fill
            subplot(3,3,7);
            plot( 1:length(chosen_top_fills), chosen_top_fills, '-o' );
            title('chosen datapoints, top fill');
            makeFigureNice( OPTS );
            % chosen bottom fill
            subplot(3,3,8);
            plot( 1:length(chosen_bot_fills), chosen_bot_fills, '-o' );
            title('chosen datapoints, bottom fill');
            makeFigureNice( OPTS );
            % chosen directivity
            subplot(3,3,9);
            plot( 1:length(chosen_dirs), 10*log10(chosen_dirs), '-o' );
            title('chosen datapoints, directivity (dB)');
            makeFigureNice( OPTS );
            
            % apodization starts here
            
            % find the optimal starting point on the power vs. pos curve
            [ obj, best_start_power ] = obj.optimize_start_power( desired_field, chosen_params );
            
            % grab and save final datapoints
            [ obj, synthesized_design ]         = obj.pick_final_datapoints_v2(desired_field, best_start_power, chosen_params);
            obj.synthesized_design              = catstruct( obj.synthesized_design, synthesized_design );  % in aux function
            obj.synthesized_design.input_wg_type = input_wg_type; 
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end function generate_final_design_apodized_v2_20200612()
        
        
        function obj = generate_final_design_apodized_gaussian( obj, MFD, input_wg_type, ...
                enforce_min_feat_size_func, fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % Synthesizes an apodized grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'bottom' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill )   
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, apodized
            if nargin < 5
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func );
            elseif nargin < 7
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override );
            else
                obj = obj.generate_final_design_apodized( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override, fill_vs_top_bot_both );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
            
        end % end function generate_final_design_apodized_gaussian()
        
        function obj = generate_final_design_apodized_gaussian_v2_20200612( obj, MFD, input_wg_type, ...
                enforce_min_feat_size_func, fill_top_override, fill_bot_override, fill_vs_top_bot_both )
            % Synthesizes an apodized grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'bottom' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill )   
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL: fill threshold override
            %   fill_vs_top_bot_both
            %       type: string
            %       desc: OPTIONAL: pick design vs top fill, bottom fill,
            %               or both
            %             either 'top', 'bottom', or 'both'
            %             default to 'top'
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, apodized
            if nargin < 5
                obj = obj.generate_final_design_apodized_v2_20200612( field_profile, input_wg_type, enforce_min_feat_size_func );
            elseif nargin < 7
                obj = obj.generate_final_design_apodized_v2_20200612( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override );
            else
                obj = obj.generate_final_design_apodized_v2_20200612( field_profile, input_wg_type, enforce_min_feat_size_func, ...
                                                            fill_top_override, fill_bot_override, fill_vs_top_bot_both );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
            
        end % end function generate_final_design_apodized_gaussian_v2_20200612()
        
        function obj = generate_final_design_uniform( obj, desired_field, input_wg_type, enforce_min_feat_size_func, filltop, fillbot )
            % Synthesizes a uniform grating that has amplitude matching field_profile
            %
            % reworking this code?
            % inputs:
            %   desired_field
            %       type: double, array
            %       desc: field profile to amplitude match to. dimensions
            %             are amplitude vs. x coordinates
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'none' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL, user defined function for enforcing minimum feature sizes
            %               currently this function must take 4 args - ( period, top_fill, bot_fill, offset )
            %   filltop
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            %   fillbot
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            
                    
            % pick the best fills
            % select subset of domain to use
            % dimensions are bot vs. top
            if nargin < 4
                [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks ] ...
                           = obj.pick_design_quadrant( input_wg_type );
            else
                % override fills
                [ obj, top_fills, bot_fills, directivities, ...
                           angles, periods, offsets, scatter_strs, ks ] ...
                           = obj.pick_design_quadrant( input_wg_type , filltop, fillbot );
            end
                   
            
            % init saving vars
            % dims are bot fill vs top fill
            obj.sweep_variables.overlap_predict_vs_fill = zeros( size( top_fills ) );

            % calculate transmission, bot fill vs top fill
            t = directivities./(directivities+1);

            % for each fill
            for i_bot = 1:size( bot_fills, 1 )

                for i_top = 1:size( top_fills, 2 )

                    % grab current alpha
                    cur_alpha = scatter_strs( i_bot, i_top );

                    % calc length of grating to scatter 99% of light (1/e^4 power)
                    grat_len = 2/cur_alpha;     % in 'units'

                    % gen x coords
                    xvec = 0 : obj.discretization : grat_len - obj.discretization;

                    % integrate the picked scatter strengths to get the predicted
                    % field shape
                    scatter_str_vs_x                = cur_alpha .* ones( size(xvec) );
                    [obj, field_shape_prediction]   = obj.predict_field_shape( xvec, scatter_str_vs_x );

                    % calculate overlap with desired field (simple xcorr)
                    [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
                    obj.sweep_variables.overlap_predict_vs_fill(i_bot, i_top) = max( abs(overlap) );

                    % check for min feat size, if doesn't satisfy, set overlap to 0
                    if exist( 'enforce_min_feat_size_func', 'var' )
                        if ~enforce_min_feat_size_func( periods( i_bot, i_top ), ...
                                top_fills( i_bot, i_top ), ...
                                bot_fills( i_bot, i_top ), ...
                                offsets( i_bot, i_top ) )
                            obj.sweep_variables.overlap_predict_vs_fill(i_bot, i_top) = 0;
                        end
                    end

                end     % end for i_top

            end     % end for i_bot

            % calculate coupling
            obj.sweep_variables.coupling_predict_vs_fill = obj.sweep_variables.overlap_predict_vs_fill.*t;
            obj.sweep_variables.fill_top_predict = top_fills( 1, : );
            obj.sweep_variables.fill_bot_predict = bot_fills( :, 1 );

            % now pick the best datapoint
            [ ~, i_bot_best, i_top_best ] = f_get_max_matrix_elem( obj.sweep_variables.coupling_predict_vs_fill );

            % generate final design
            grat_len            = 4./scatter_strs( i_bot_best, i_top_best );
            xvec                = 0 : obj.discretization : grat_len - obj.discretization;
            alpha_des           = scatter_strs( i_bot_best, i_top_best ) .* ones( size( xvec ) );
            [ obj, synthesized_des ] = obj.pick_final_datapoints( xvec, alpha_des, ...
                                             directivities( i_bot_best, i_top_best ), ...
                                             bot_fills( i_bot_best, i_top_best ), ...
                                             top_fills( i_bot_best, i_top_best ), ...
                                             offsets( i_bot_best, i_top_best ), ...
                                             periods( i_bot_best, i_top_best ), ...
                                             angles( i_bot_best, i_top_best ), ...
                                             scatter_strs( i_bot_best, i_top_best ), ...
                                             ks( i_bot_best, i_top_best ) );
            obj.synthesized_design = synthesized_des;
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % build final index distribution
            obj = obj.build_final_index();
            
        end     % end generate_final_design_uniform()
        
        
                    
        function obj = generate_final_design_uniform_gaussian( obj, MFD, input_wg_type, enforce_min_feat_size_func, filltop, fillbot )
            % Synthesizes a uniform grating that radiates desired Gaussian field profile
            %
            % Inputs:
            %   MFD
            %       type: double, scalar
            %       desc: mode field diameter of gaussian, defined as 1/e^2 width of intensity
            %   input_wg_type
            %       type: string
            %       desc: 'full' or 'none' are currently supported
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: user defined function for enforcing minimum feature sizes
            %               currently this function must take 2 args - ( period, fill ) 
            %   filltop
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            %   fillbot
            %       type: double, scalar
            %       desc: OPTIONAL fill to use
            
            % generate a fiber gaussian mode
            [ obj, field_profile ] = obj.make_gaussian_profile( MFD );
            
            % generate final design, uniform
            if nargin > 4
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type, enforce_min_feat_size_func, filltop, fillbot );
            else
                obj = obj.generate_final_design_uniform( field_profile, input_wg_type, enforce_min_feat_size_func );
            end
            
            % save MFD
            obj.synthesized_design.MFD = MFD;
        
        end     % end generate_final_design_uniform_gaussian()
        
        function obj = choose_unit_cells( obj, input_wg_type, enforce_min_feat_size_func, ...
                                        fill_top_override, fill_bot_override, fill_vs_top_bot_both, ...
                                        interp_on, plot_on )
            % added a new input - interp_on. optionally interpolates the
            % final chosen cells. default is true
            % plot_on - optional. decide whether to plot the chosen cells
            % or not. default is true
                                    
            % select subset of domain to use
            if nargin < 4
             [ obj, top_fills, bot_fills, directivities, ...
                       angles, periods, offsets, scatter_strs, ks ] ...
                       = obj.pick_design_quadrant( input_wg_type );
            else
             % override the fills as user desired
             [ obj, top_fills, bot_fills, directivities, ...
                       angles, periods, offsets, scatter_strs, ks ] ...
                       = obj.pick_design_quadrant( input_wg_type, fill_top_override, fill_bot_override );
            end
            
            if nargin < 7
                interp_on = true;
            end
            
            if nargin < 8
                plot_on = true;
            end

%             % try interpolating
%             n_interp            = 200;
%             top_fills_interp 	= linspace(min(top_fills(:)), max(top_fills(:)), n_interp);
%             bot_fills_interp    = linspace(min(bot_fills(:)), max(bot_fills(:)), n_interp);
%             [top_fills_interp, bot_fills_interp] = meshgrid(top_fills_interp, bot_fills_interp);
%             directivities       = interp2( top_fills, bot_fills, directivities, top_fills_interp, bot_fills_interp );
%             angles              = interp2( top_fills, bot_fills, angles, top_fills_interp, bot_fills_interp );
%             periods             = interp2( top_fills, bot_fills, periods, top_fills_interp, bot_fills_interp );
%             offsets             = interp2( top_fills, bot_fills, offsets, top_fills_interp, bot_fills_interp );
%             scatter_strs        = interp2( top_fills, bot_fills, scatter_strs, top_fills_interp, bot_fills_interp );
%             ks                  = interp2( top_fills, bot_fills, ks, top_fills_interp, bot_fills_interp );
%             top_fills           = top_fills_interp;
%             bot_fills           = bot_fills_interp;
            
            % default to using top fill
            if nargin < 6
                fill_vs_top_bot_both = 'top';
            end

            % for each top fill, pick bottom fill with highest
            % directivity
            [ chosen_dirs_vs_top, indx_highest_dir_per_topfill ] = max( directivities, [], 1 );
            % gotta use linear indexing
            indxs_vs_top               = sub2ind( size(directivities), indx_highest_dir_per_topfill, 1:size(directivities,2)  );


            % for each bottom fill, pick top fill
            % with highest directivity
            % directivity
            [ chosen_dirs_vs_bot, indx_highest_dir_per_botfill ] = max( directivities, [], 2 );
            % gotta use linear indexing
            indxs_vs_bot                   = sub2ind( size(directivities), 1:size(directivities,1), indx_highest_dir_per_botfill.' );

            switch fill_vs_top_bot_both

                case 'top'

                    % use top fill
                    chosen_angles       = angles( indxs_vs_top );
                    chosen_periods      = periods( indxs_vs_top );
                    chosen_offsets      = offsets( indxs_vs_top );
                    chosen_scatter_str  = scatter_strs( indxs_vs_top );
                    chosen_ks           = ks( indxs_vs_top );
                    chosen_top_fills    = top_fills( indxs_vs_top );
                    chosen_bot_fills    = bot_fills( indxs_vs_top );
                    chosen_dirs         = chosen_dirs_vs_top;

                case 'bottom'

                    % use bot fill
                    chosen_angles       = angles( indxs_vs_bot );
                    chosen_periods      = periods( indxs_vs_bot );
                    chosen_offsets      = offsets( indxs_vs_bot);
                    chosen_scatter_str  = scatter_strs( indxs_vs_bot );
                    chosen_ks           = ks( indxs_vs_bot );
                    chosen_top_fills    = top_fills( indxs_vs_bot );
                    chosen_bot_fills    = bot_fills( indxs_vs_bot );
                    chosen_dirs         = chosen_dirs_vs_bot;

                case 'both'

                    % append bot fill
                    chosen_angles       = [ angles( indxs_vs_top ), angles( indxs_vs_bot ) ];
                    chosen_periods      = [ periods( indxs_vs_top ), periods( indxs_vs_bot ) ];
                    chosen_offsets      = [ offsets( indxs_vs_top ), offsets( indxs_vs_bot) ];
                    chosen_scatter_str  = [ scatter_strs( indxs_vs_top ), scatter_strs( indxs_vs_bot ) ];
                    chosen_ks           = [ ks( indxs_vs_top ), ks( indxs_vs_bot ) ];
                    chosen_top_fills    = [ top_fills( indxs_vs_top ), top_fills( indxs_vs_bot ) ];
                    chosen_bot_fills    = [ bot_fills( indxs_vs_top ), bot_fills( indxs_vs_bot ) ];
                    chosen_dirs = [ chosen_dirs_vs_top, chosen_dirs_vs_bot.' ];

            end

            % try sorting on scattering strength
            [ chosen_scatter_str, indx_sort ] = sort( chosen_scatter_str );
            chosen_angles       = chosen_angles( indx_sort );
            chosen_periods      = chosen_periods( indx_sort );
            chosen_offsets      = chosen_offsets( indx_sort );
            chosen_ks           = chosen_ks( indx_sort );
            chosen_top_fills    = chosen_top_fills( indx_sort );
            chosen_bot_fills    = chosen_bot_fills( indx_sort );
            chosen_dirs         = chosen_dirs( indx_sort );
            
            % remove 0 directionality cells
            chosen_scatter_str  = chosen_scatter_str( chosen_dirs > 0 );
            chosen_angles       = chosen_angles( chosen_dirs > 0 );
            chosen_periods      = chosen_periods( chosen_dirs > 0 );
            chosen_offsets      = chosen_offsets( chosen_dirs > 0 );
            chosen_ks           = chosen_ks( chosen_dirs > 0 );
            chosen_top_fills    = chosen_top_fills( chosen_dirs > 0 );
            chosen_bot_fills    = chosen_bot_fills( chosen_dirs > 0 );
            chosen_dirs         = chosen_dirs( chosen_dirs > 0 );
            
            % how about interpolating here
            if interp_on
                n_interp            = 200;
                chosen_scatter_str  = interp1( 1:length(chosen_scatter_str), chosen_scatter_str, linspace(1, length(chosen_scatter_str), n_interp) );
                chosen_angles       = interp1( 1:length(chosen_angles), chosen_angles, linspace(1, length(chosen_angles), n_interp) );
                chosen_periods      = interp1( 1:length(chosen_periods), chosen_periods, linspace(1, length(chosen_periods), n_interp) );
                chosen_offsets      = interp1( 1:length(chosen_offsets), chosen_offsets, linspace(1, length(chosen_offsets), n_interp) );
                chosen_ks           = interp1( 1:length(chosen_ks), chosen_ks, linspace(1, length(chosen_ks), n_interp) );
                chosen_top_fills    = interp1( 1:length(chosen_top_fills), chosen_top_fills, linspace(1, length(chosen_top_fills), n_interp) );
                chosen_bot_fills    = interp1( 1:length(chosen_bot_fills), chosen_bot_fills, linspace(1, length(chosen_bot_fills), n_interp) );
                chosen_dirs         = interp1( 1:length(chosen_dirs), chosen_dirs, linspace(1, length(chosen_dirs), n_interp) );
            end
                
            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;

                indices_to_keep = [];
                for ii = 1:length( chosen_periods )
                    if enforce_min_feat_size_func( chosen_periods(ii), chosen_top_fills(ii), chosen_bot_fills(ii), chosen_offsets(ii) ) == true
                        indices_to_keep(end+1) = ii;
                    end
                end

                chosen_angles       = chosen_angles(indices_to_keep);
                chosen_periods      = chosen_periods(indices_to_keep);
                chosen_offsets      = chosen_offsets(indices_to_keep);
                chosen_scatter_str  = chosen_scatter_str(indices_to_keep);
                chosen_ks           = chosen_ks(indices_to_keep);
                chosen_top_fills    = chosen_top_fills(indices_to_keep);
                chosen_bot_fills    = chosen_bot_fills(indices_to_keep);
                chosen_dirs         = chosen_dirs(indices_to_keep);

            end
           
            if plot_on
                OPTS = struct('fig_size', [2400, 1200]);
                figure('Name', 'chosen_datapoints');
                % chosen angles
                subplot(3,3,1);
                plot( 1:length(chosen_angles), chosen_angles, '-o' );
                title('chosen datapoints, angles');
                xlim( [1, length(chosen_angles) ] );
                makeFigureNice( OPTS );
                % chosen periods
                subplot(3,3,2);
                plot( 1:length(chosen_periods), chosen_periods, '-o' );
                title('chosen datapoints, periods');
                makeFigureNice( OPTS );
                % chosen offsets
                subplot(3,3,3);
                plot( 1:length(chosen_offsets), chosen_offsets, '-o' );
                title('chosen datapoints, offsets');
                makeFigureNice( OPTS );
                % chosen scatter strength
                subplot(3,3,4);
                plot( 1:length(chosen_scatter_str), chosen_scatter_str, '-o' );
                title('chosen datapoints, scattering strengths');
                makeFigureNice( OPTS );
                % chosen k real
                subplot(3,3,5);
                plot( 1:length(chosen_ks), real(chosen_ks), '-o' );
                title('chosen datapoints, k real');
                makeFigureNice( OPTS );
                % chosen k imag
                subplot(3,3,6);
                plot( 1:length(chosen_ks), imag(chosen_ks), '-o' );
                title('chosen datapoints, k imaginary');
                makeFigureNice( OPTS );
                % chosen top fill
                subplot(3,3,7);
                plot( 1:length(chosen_top_fills), chosen_top_fills, '-o' );
                title('chosen datapoints, top fill');
                makeFigureNice( OPTS );
                % chosen bottom fill
                subplot(3,3,8);
                plot( 1:length(chosen_bot_fills), chosen_bot_fills, '-o' );
                title('chosen datapoints, bottom fill');
                makeFigureNice( OPTS );
                % chosen directivity
                subplot(3,3,9);
                plot( 1:length(chosen_dirs), 10*log10(chosen_dirs), '-o' );
                title('chosen datapoints, directivity (dB)');
                makeFigureNice( OPTS );
            end

            % save the chosen cells
            obj.chosen_cells = struct( 'chosen_angles', chosen_angles, ...
                                    'chosen_periods', chosen_periods, ...
                                    'chosen_offsets', chosen_offsets, ...
                                    'chosen_scatter_str', chosen_scatter_str, ...
                                    'chosen_ks', chosen_ks, ...
                                    'chosen_top_fills', chosen_top_fills, ...
                                    'chosen_bot_fills', chosen_bot_fills, ...
                                    'chosen_dirs', chosen_dirs );
                                    
        end % end choose_unit_cells()
        
       
        
        function [ obj, best_alpha_power ] = optimize_start_alpha( obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field )
            % Optimizes the starting alpha, based on predicted field overlap
            %
            % Inputs:
            %     desired_field
            %         type: double, array
            %         desc: desired field vs. xvec
            
%             % pick the alphas to try:
%             alpha_power_max = log10( min( scatter_strs_high_dir ) );
%             alpha_powers    = linspace( -6, alpha_power_max, 100 );
%             alpha_try       = 10.^( alpha_powers );
            
%             % save max overlap, dimensions overlap vs. alpha
%             max_overlap = zeros( size( alpha_try ) );
            
%             % brute force optimize
%             for i_alpha = 1:length(alpha_powers)
%                 
%                 max_overlap(i_alpha) = obj.predict_overlap_for_optimization( ...
%                                               xvec, alpha_des, high_dirs, ...
%                                               bot_fills_high_dir, top_fills_high_dir, ...
%                                               offsets_high_dir, periods_high_dir, ...
%                                               angles_high_dir, scatter_strs_high_dir, ...
%                                               k_high_dir, desired_field, alpha_powers(i_alpha) );
%                 
%             end     % end for i_alpha = 1:length(alpha_try)
           
            % using matlab's fminsearch
            f = @(alpha_power)( 1 - obj.predict_overlap_for_optimization( ...
                                              xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field, alpha_power ) );
            
            % options
%             opts = optimset( 'Display', 'iter' );
                                          
            % for debuggging
%             [ best_alpha_power, best_val] = fminsearch( f, log10(min( scatter_strs_high_dir )), opts );
            
            best_alpha_power = fminsearch( f, log10(min( chosen_scatter_strs )) );
            
%             best_overlap_fminsearch = 1 - best_val;
            
%             % DEBUG plot max overlap vs. alpha
%             figure;
%             plot( alpha_powers, max_overlap, '-o' );
%             xlabel('alpha power'); ylabel('max overlap');
%             title('DEBUG max overlap vs. start alpha (in power)');
%             makeFigureNice();
            
            % now save the final design with the best predicted overlap
            obj = obj.pick_final_datapoints(  xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, 10.^(best_alpha_power) );
                                          
        end     % end function optimize_start_alpha()
        
        function [ obj, best_start_power ] = optimize_start_power( obj, desired_field, chosenparams )
            % Optimizes the starting radiation power, based on predicted field overlap
            % for the new version of apodization, 6/14/2020
            %
            % Inputs:
            %     desired_field
            %         type: double, array
            %         desc: desired field vs. xvec
           
            % using matlab's fminsearch
            f = @(startpower)( 1 - obj.predict_overlap_for_optimization_v2( desired_field, startpower, chosenparams ));                            
            best_start_power = fminsearch( f, min(chosenparams.chosen_rad_power_ratio) );
                                                     
        end     % end function optimize_start_alpha()
        
        
        function [ max_overlap, field_shape_prediction ] = predict_overlap_for_optimization( ...
                                              obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, desired_field, start_alpha_power )
            % merit function used in optimize_start_alpha for predicting
            % overlap
            %
            % Inputs:
            %   start_alpha_power:
            %       type: double, scalar
            %       desc: alpha to try = 10^(start_alpha_power)
            %
            % Outputs:
            %     max_overlap
            %     field_shape_prediction
            %         type: double, array
            %         desc: predicted field shape vs. x, mostly for debugging
                      
            % pick the datapoints
            [ obj, synthesized_design ] = obj.pick_final_datapoints(  xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, 10.^(start_alpha_power) );
                                          
            % convert scatter str vs. cell into scatter str vs. x
            xvec                = xvec - xvec(1);
            scatter_str_vs_x    = zeros( size( xvec ) );
            end_pos_vs_cell     = cumsum( synthesized_design.period );
            cur_cell            = 1;
            for ii = 1:length(xvec)

                scatter_str_vs_x(ii) = synthesized_design.scatter_str(cur_cell);

                if xvec(ii) > end_pos_vs_cell( cur_cell ) && cur_cell < length(end_pos_vs_cell)
                    % move to next cell
                    cur_cell = cur_cell + 1;
                end

            end

            % integrate the picked scatter strengths to get the predicted
            % field shape
            
            [obj, field_shape_prediction] = obj.predict_field_shape( xvec, scatter_str_vs_x );
            
            % calculate overlap with desired field
            [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
            max_overlap     = max( abs(overlap) );
            
        end     % end function predict_overlap_for_optimization()
  
        function [ max_overlap, field_shape_prediction ] = predict_overlap_for_optimization_v2( ...
                                              obj, desired_field, starting_power, chosenparams )
            % merit function used in optimize_start_alpha for predicting
            % overlap   
           
            [ obj, synthesized_design ] = obj.pick_final_datapoints_v2(desired_field, starting_power, chosenparams);
                                          
            % convert power ratio into scatter strength vs cell
            scatter_str = (-1./(2*synthesized_design.period)) .* log( 1 - synthesized_design.chosen_rad_power_ratio );
            
            % convert scatter str vs. cell into scatter str vs. x
            xvec                = synthesized_design.xvec - synthesized_design.xvec(1);
            scatter_str_vs_x    = zeros( size( xvec ) );
            imag_k_vs_x         = zeros( size( xvec ) );
            end_pos_vs_cell     = cumsum( synthesized_design.period );
            cur_cell            = 1;
            for ii = 1:length(xvec)

                scatter_str_vs_x(ii) = scatter_str(cur_cell);
                imag_k_vs_x(ii)     = imag(synthesized_design.k(cur_cell));

                if xvec(ii) > end_pos_vs_cell( cur_cell ) && cur_cell < length(end_pos_vs_cell)
                    % move to next cell
                    cur_cell = cur_cell + 1;
                end

            end

            % integrate the picked scatter strengths to get the predicted
            % field shape
            
            [obj, field_shape_prediction] = obj.predict_field_shape_v2( xvec, scatter_str_vs_x, imag_k_vs_x );
            
            % calculate overlap with desired field
            [ ~, overlap ]  = obj.xcorr_normalized( desired_field, field_shape_prediction );
            max_overlap     = max( abs(overlap) );
            
        end     % end function predict_overlap_for_optimization_v2()
  
        
        
        function [ obj, synthesized_design ] = pick_final_datapoints( obj, xvec, alpha_des, chosen_dirs, ...
                                              chosen_bot_fills, chosen_top_fills, ...
                                              chosen_offsets, chosen_periods, ...
                                              chosen_angles, chosen_scatter_strs, ...
                                              chosen_ks, start_alpha_des )
            % Picks the final datapoints (cells) that make up the grating
            %
            % Inputs:
            %     xvec
            %         type: double, array
            %         desc: coordinates of desired field profile
            %     alpha_des
            %         type: double, array
            %         desc: desired scattering strength vs. x
            %     high_dirs
            %         type: double, array
            %         desc: chosen highest directivity points
            %     bot_fills_high_dir
            %         type: double, array
            %         desc: bottom fill value of chosen highest dir. points
            %     top_fills_high_dir
            %         type: double, array
            %         desc: top fill value of chosen highest dir. points
            %     offsets_high_dir
            %         type: double, array
            %         desc: offset values of chosen highest dir. points
            %     periods_high_dir
            %         type: double, array
            %         desc: periods of chosen highest dir. points
            %     angles_high_dir
            %         type: double, array
            %         desc: angle of chosen highest dir. points
            %     scatter_strs_high_dir
            %         type: double, array
            %         desc: scattering strength (alpha) of chosen highest
            %               dir. points
            %     k_high_dir
            %         type: double, array
            %         desc: k of chosen highest dir. points
            %     start_alpha_des
            %         type: double, scalar
            %         desc: OPTIONAL, user defined starting alpha
            %
            % Outputs:
            %   synthesized_design:
            %       type: struct
            %       desc: a struct that holds the synthesized cells
            
            % default to weakest cell
            if ~exist( 'start_alpha_des', 'var' )
                start_alpha_des = min(chosen_scatter_strs); %1e-5;
            end
            
            % now match these data points to the desired alpha
            % starting point
            [~, indx_max_alpha] = max( alpha_des );
            [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
            cur_x               = xvec(indx_x);
            
            % final synthesized variables
            synthesized_design.dir                  = [];
            synthesized_design.bot_fill             = [];
            synthesized_design.top_bot_fill_ratio   = [];
            synthesized_design.top_fill             = [];
            synthesized_design.period               = [];
            synthesized_design.offset               = [];
            synthesized_design.angles               = [];
            synthesized_design.scatter_str          = [];
            synthesized_design.k                    = [];
            synthesized_design.GC                   = {};
            synthesized_design.des_scatter          = [];
            
            
            % flag for switching to using max scattering strength
            saturate_scatter_str_to_max = false;
 
            ii = 1;
            while cur_x < xvec(end)
                % build grating one cell at a time
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter = alpha_des(indx_x);                            % desired alpha
                if ii > 1
                    if des_scatter  > max( chosen_scatter_strs ) || des_scatter < synthesized_design.des_scatter(ii-1)
                        % desired scattering strength too high, gotta saturate
                        saturate_scatter_str_to_max = true;
                    end
                end
                if ~saturate_scatter_str_to_max
                    [~, indx_closest_scatter]   = min( abs(chosen_scatter_strs - des_scatter) );          % index of closest scatter design 
                else
                    [~, indx_closest_scatter]   = max( chosen_scatter_strs );                             % saturate to max
                end
                
                % save parameters
                synthesized_design.dir(ii)                   = chosen_dirs( indx_closest_scatter );
                synthesized_design.bot_fill(ii)              = chosen_bot_fills( indx_closest_scatter );
                synthesized_design.top_fill(ii)              = chosen_top_fills( indx_closest_scatter );
                synthesized_design.top_bot_fill_ratio(ii)    = chosen_top_fills( indx_closest_scatter ) ./ chosen_bot_fills( indx_closest_scatter );
                synthesized_design.offset(ii)                = chosen_offsets( indx_closest_scatter );
                synthesized_design.period(ii)                = chosen_periods( indx_closest_scatter );
                synthesized_design.angles(ii)                = chosen_angles( indx_closest_scatter );
                synthesized_design.scatter_str(ii)           = chosen_scatter_strs( indx_closest_scatter );
                synthesized_design.k(ii)                     = chosen_ks( indx_closest_scatter );
                synthesized_design.des_scatter(ii)           = des_scatter;
                
                synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
                               obj.discretization, ...
                               obj.background_index, ...
                               obj.y_domain_size, ...
                               synthesized_design.period(ii), ...
                               synthesized_design.top_fill(ii), ...
                               synthesized_design.bot_fill(ii), ...
                               synthesized_design.offset(ii)/synthesized_design.period(ii) );
                           
                % move onto next
                cur_x       = cur_x + synthesized_design.period(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                ii          = ii + 1;
                
            end     % end for ii = 1:ncells
            
        end     % end pick_final_datapoints()
              
        function [ obj, synthesized_design ] = pick_final_datapoints_v2(obj, desired_field, starting_power, chosenparams)
            % new version, 6/14/2020
            % follows total radiated power profile of the desired field to
            % synthesize an apodized grating
            
            % generate x coordinates for the field profile
            xvec = obj.discretization .* ( 0:length(desired_field)-1 );
            xvec = xvec - xvec(round(end/2));                                % shift origin over to middle

            % calculate desired total radiated power vs position
            power_vs_pos    = cumsum( abs(desired_field).^2 ) * obj.discretization;
            power_vs_pos    = power_vs_pos./max(abs(power_vs_pos));     % normalize power to 1
            
            % for now, pick a starting power and apodize from there
            % starting point
%             starting_power      = min(chosen_rad_power_ratio);
            [~, indx_startx]    = min(abs( power_vs_pos - starting_power ) );
            cur_x               = xvec(indx_startx);
            cur_x_vs_cell       = []; % for debugging
                      
            % pick datapoints
            
            % final synthesized variables
            synthesized_design.dir                  = [];
            synthesized_design.bot_fill             = [];
            synthesized_design.top_bot_fill_ratio   = [];
            synthesized_design.top_fill             = [];
            synthesized_design.period               = [];
            synthesized_design.offset               = [];
            synthesized_design.angles               = [];
            synthesized_design.scatter_str          = [];
            synthesized_design.k                    = [];
            synthesized_design.GC                   = {};
%             synthesized_design.des_rad_power        = [];   % for debugging
            synthesized_design.chosen_rad_power_ratio     = [];
            synthesized_design.xvec                 = xvec;
            synthesized_design.power_vs_pos         = power_vs_pos;
            
            % starting input power
            Pin = max(power_vs_pos);
            Pin_vs_cell                 = []; %  mostly for debugging
            power_rad_err_vs_cell       = []; % for debugging
            desired_rad_power_vs_cell   = [];
            pred_rad_power_vs_cell      = [];
            cur_rad_power_vs_cell       = [];
            
            
            ii = 1;
            cur_rad_power = 0;
            while cur_x < xvec(end) && cur_rad_power < 1
                % build grating one cell at a time
                
                % first go through each possible cell and find the one that
                % scatters the best
                power_rad_error = zeros(size(chosenparams.chosen_rad_power_ratio));
                pred_rad_power  = zeros(size(chosenparams.chosen_rad_power_ratio));
                desired_rad_power = zeros(size(chosenparams.chosen_rad_power_ratio));
                for i_cell = 1:length(chosenparams.chosen_periods)
                    
                    this_period     = chosenparams.chosen_periods(i_cell);
                    next_x          = cur_x + this_period;
                    [~, indx_next_power]    = min(abs(xvec - next_x));
%                     [~, indx_cur_power]     = min(abs(xvec - cur_x));
                    desired_rad_power(i_cell) = power_vs_pos( indx_next_power ) - cur_rad_power;
                    pred_rad_power(i_cell)  = chosenparams.chosen_rad_power_ratio(i_cell) * Pin;
                    power_rad_error(i_cell) = abs( pred_rad_power(i_cell) - desired_rad_power(i_cell) );
                    
                end
                
                % pick cell with radiated power closest to desired
                [power_rad_err_vs_cell(ii), indx_closest_scatter ] = min(power_rad_error);
                desired_rad_power_vs_cell(ii)   = desired_rad_power(indx_closest_scatter);
                pred_rad_power_vs_cell(ii)      = pred_rad_power(indx_closest_scatter);
                
                % update current amount of radiated power
                cur_rad_power_vs_cell(ii) = cur_rad_power;
                cur_rad_power = cur_rad_power + pred_rad_power_vs_cell(ii);
                
                
                % save parameters
                synthesized_design.dir(ii)                   = chosenparams.chosen_dirs( indx_closest_scatter );
                synthesized_design.bot_fill(ii)              = chosenparams.chosen_bot_fills( indx_closest_scatter );
                synthesized_design.top_fill(ii)              = chosenparams.chosen_top_fills( indx_closest_scatter );
                synthesized_design.top_bot_fill_ratio(ii)    = chosenparams.chosen_top_fills( indx_closest_scatter ) ./ chosenparams.chosen_bot_fills( indx_closest_scatter );
                synthesized_design.offset(ii)                = chosenparams.chosen_offsets( indx_closest_scatter );
                synthesized_design.period(ii)                = chosenparams.chosen_periods( indx_closest_scatter );
                synthesized_design.angles(ii)                = chosenparams.chosen_angles( indx_closest_scatter );
                synthesized_design.scatter_str(ii)           = chosenparams.chosen_scatter_str( indx_closest_scatter );
                synthesized_design.k(ii)                     = chosenparams.chosen_ks( indx_closest_scatter );
                synthesized_design.chosen_rad_power_ratio(ii)      = chosenparams.chosen_rad_power_ratio( indx_closest_scatter );
%                 synthesized_design.des_rad_power(ii)         = des_scatter;
                
                synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
                               obj.discretization, ...
                               obj.units.name, ...
                               obj.lambda, ...
                               obj.background_index, ...
                               obj.y_domain_size, ...
                               synthesized_design.period(ii), ...
                               synthesized_design.top_fill(ii), ...
                               synthesized_design.bot_fill(ii), ...
                               synthesized_design.offset(ii)/synthesized_design.period(ii) );
                           
                % move onto next
                cur_x_vs_cell(ii) = cur_x;
                cur_x       = cur_x + synthesized_design.period(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                Pin_vs_cell(ii) = Pin;
%                 Pin         = Pin * exp( -2*imag(synthesized_design.k(ii))*synthesized_design.period(ii) );
                Pin         = Pin - pred_rad_power_vs_cell(ii);
                ii          = ii + 1;
                
            end     % end cur_x < xvec(end)
            
            synthesized_design.cur_x_vs_cell = cur_x_vs_cell;
            synthesized_design.Pin_vs_cell = Pin_vs_cell;
            synthesized_design.power_rad_err_vs_cell = power_rad_err_vs_cell;
            synthesized_design.desired_rad_power_vs_cell = desired_rad_power_vs_cell;
            synthesized_design.pred_rad_power_vs_cell = pred_rad_power_vs_cell;
            synthesized_design.pred_total_rad_power_vs_cell = cumsum(pred_rad_power_vs_cell);
            synthesized_design.des_total_rad_power_vs_cell = cumsum(desired_rad_power_vs_cell);
            synthesized_design.cur_rad_power_vs_cell = cur_rad_power_vs_cell;
            
        end     % end pick_final_datapoints_v2
        
        function [ obj, top_fills, bot_fills, directivities, ...
                   angles, periods, offsets, scatter_strs, ks, rad_power_ratio ] ...
                   = pick_design_quadrant( obj, input_wg_type, fill_top_override, fill_bot_override )
            % Picks out the design quadrant of the top/bot fill space that
            % corresponds to the desired input waveguide type
            %
            % Inputs:
            %   input_wg_type
            %       type: string
            %       desc: 'top' 'bottom' 'full' 'none'
            %   fill_top_override
            %       type: scalar
            %       desc: OPTIONAL - overrides the center of mass
            %   fill_bot_override
            %       type: scalar
            %       desc: OPTIONAL - overrides the center of mass 
            
            % meshgrid the fills (dimensions are bot vs. top, I believe)
            [ top_fills_mesh, bot_fills_mesh ] = meshgrid( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots );

            % find where the max scattering strength is
            if nargin < 3
                % new ver. where i find center of mass
                r       = [ top_fills_mesh(:), bot_fills_mesh(:) ];   % dimensions [ top fill vs. mass, bot fill vs. mass ]
                mass    = obj.sweep_variables.scatter_str_vs_fills(:); 
                R       = sum( r .* [ mass, mass ], 1 )./sum(mass(:));
                [~, i_max_top_scat] = min(abs( obj.sweep_variables.fill_tops - R(1) ) );
                [~, i_max_bot_scat] = min(abs( obj.sweep_variables.fill_bots - R(2) ) );
            else
                % override center of mass
                [~, i_max_top_scat] = min(abs( obj.sweep_variables.fill_tops - fill_top_override ) );
                [~, i_max_bot_scat] = min(abs( obj.sweep_variables.fill_bots - fill_bot_override ) );
            end
            
            % Select subset of domain to use
            if strcmp( input_wg_type, 'bottom' ) == true
                % input waveguide is the bottom layer
            
                % narrow down design space (take top left quadrant)
                % recall fills are sorted in descending order
                top_fills       = top_fills_mesh( 1:i_max_bot_scat, i_max_top_scat:end );
                bot_fills       = bot_fills_mesh( 1:i_max_bot_scat, i_max_top_scat:end );
                directivities   = obj.sweep_variables.directivities_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                angles          = obj.sweep_variables.angles_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );      
                periods         = obj.sweep_variables.periods_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                offsets         = obj.sweep_variables.offsets_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                ks              = obj.sweep_variables.k_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( 1:i_max_bot_scat, i_max_top_scat:end );
                catch
                end
                
            elseif strcmp( input_wg_type, 'full' ) == true
                % input waveguide is a full waveguide
                
                % narrow down design space (take top right quadrant)
                top_fills       = top_fills_mesh( 1:i_max_bot_scat, 1:i_max_top_scat );
                bot_fills       = bot_fills_mesh( 1:i_max_bot_scat, 1:i_max_top_scat );
                directivities   = obj.sweep_variables.directivities_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                angles          = obj.sweep_variables.angles_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );      
                periods         = obj.sweep_variables.periods_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                offsets         = obj.sweep_variables.offsets_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                ks              = obj.sweep_variables.k_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( 1:i_max_bot_scat, 1:i_max_top_scat );
                catch
                end
                    
            elseif strcmp( input_wg_type, 'top' ) == true
                % input waveguide is a full waveguide
                
                % narrow down design space (take bottom right quadrant)
                top_fills       = top_fills_mesh( i_max_bot_scat:end, 1:i_max_top_scat );
                bot_fills       = bot_fills_mesh( i_max_bot_scat:end, 1:i_max_top_scat );
                directivities   = obj.sweep_variables.directivities_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                angles          = obj.sweep_variables.angles_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );      
                periods         = obj.sweep_variables.periods_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                offsets         = obj.sweep_variables.offsets_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                ks              = obj.sweep_variables.k_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( i_max_bot_scat:end, 1:i_max_top_scat );
                catch
                end
                
            elseif strcmp( input_wg_type, 'none' ) == true
                % input waveguide is cladding
                
                % narrow down design space (take bottom left quadrant)
                top_fills       = top_fills_mesh( i_max_bot_scat:end, i_max_top_scat:end );
                bot_fills       = bot_fills_mesh( i_max_bot_scat:end, i_max_top_scat:end );
                directivities   = obj.sweep_variables.directivities_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                angles          = obj.sweep_variables.angles_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );      
                periods         = obj.sweep_variables.periods_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                offsets         = obj.sweep_variables.offsets_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                ks              = obj.sweep_variables.k_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                try
                    rad_power_ratio = obj.sweep_variables.prad_pin_vs_fills( i_max_bot_scat:end, i_max_top_scat:end );
                catch
                end
                
            end     % end if strcmp( input_wg_type, 'bottom' )
            
        end     % end pick_design_quadrant()

        
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
        
        
        function obj = synthesizeGaussianGrating(obj, MFD, DEBUG)
            % DEPRECATED
            % replaced (mostly) by generate_design_space
            %
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
            figure('Name','des_scatter');
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice();
            
            % -------------------------------------------------------------
            % Simulation time
            
            
           
            % get waveguide k
            fprintf('Simulating waveguide...\n');           

            
            % ----------------------------------------
            % NEW NEW VERSION SWEEPING JELENAS DATAPOINTS
            % SPLITTING NORMAL AND INVERTED DOMAIN
            % ----------------------------------------
            
            fprintf('Sweeping fill factors for directivity and angle...\n');
            
            % set fill factors and offsets
%             fill_bots           = fliplr( 0.4:0.025:1.0 );
%             fill_top_bot_ratio  = fliplr( 0.0:0.025:1.2 );
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
















































