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
        % offset              
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
        
    end
    
    methods
        
        function obj = c_synthTwoLevelGrating(varargin)
            % Constructor
            % See top comments for input documentation
            
            % call c_synthGrating constructor
            obj = obj@c_synthGrating(varargin{:});
            
            % inputs and defaults
            inputs      = { 'coupling_direction',   'none' }; 
            obj.inputs  = [ obj.inputs, inputs ];                           % append inputs
            
            % parse inputs
            p = f_parse_varargin( inputs, varargin{:} );
            
            obj.coupling_direction = p.coupling_direction;

        end     % end constructor()

        
        function obj = generate_design_space( obj, fill_bots, fill_top_bot_ratio, verbose )
            % Generates the design space vs. fill factors
            % picking optimum offset and period for highest directivity and
            % closest angle to desired
            %
            % based on synthesizeGaussianGrating
            %
            % would be good to implement: option to save GC data or not
            %
            % Inputs:
            %   fill_bots
            %       type: double, array
            %       desc: OPTIONAL Currently mostly for testing
            %   fill_top_bot_ratio
            %       type: double, array
            %       desc: OPTIONAL Currently mostly for testing
            %   verbose
            %       type: double, array
            %       desc: OPTIONAL Currently mostly for testing, spits out
            %             a bunch of stuff to the prompt
            
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
                fill_bots           = fliplr( 0.4:0.025:0.975 );
            end
            if ~exist('fill_top_bot_ratio', 'var')
                % default fill
                fill_top_bot_ratio  = fliplr( 0.05:0.025:1.2 );
            end
            %             fill_top_bot_ratio  = fliplr( 0.9:0.025:1.2 );
%             fill_bots           = fliplr( 0.95:0.025:0.975 );
%             fill_top_bot_ratio  = fliplr( 0.95:0.025:1 );
%             fill_bots           = fliplr( 0.475:0.025:0.975 );
%             fill_top_bot_ratio  = fliplr( 1.175:0.025:1.2 );
            fill_tops           = [];                                       %fill_bots .* fill_top_bot_ratio;
            guess_offset        = 0;
            
            % sort fills so they are in descending order
            fill_bots           = sort( fill_bots, 'descend' );
            fill_top_bot_ratio  = sort( fill_top_bot_ratio, 'descend' );
            
            % save fills and offsets
            obj.sweep_variables.fill_tops           = fill_tops;
            obj.sweep_variables.fill_bots           = fill_bots;
            obj.sweep_variables.fill_top_bot_ratio  = fill_top_bot_ratio;
            
            % initialize saving variables
            directivities_vs_fills  = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            angles_vs_fills         = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            periods_vs_fills        = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            offsets_vs_fills        = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio, this is offset ratio
            scatter_str_vs_fills    = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            k_vs_fills              = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            GC_vs_fills             = cell( length( fill_bots ), length( fill_top_bot_ratio ) );      % dimensions bot fill vs. top/bot ratio
            dir_b4_period_vs_fills  = zeros( length( fill_bots ), length( fill_top_bot_ratio ) );     % dimensions bot fill vs. top/bot ratio
            
            % make grating cell, assuming both layers are filled
            waveguide = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.units.name, ...
                                               obj.lambda, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               2*obj.discretization, ...
                                               1.0, 1.0, 0.0 );
            
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
                fill_top = fill_top_bot_ratio(1) * fill_bots(i_ff_bot);
                if fill_top < 1
                    % Only run optimization if theres a perturbation of
                    % both layers

                    % run optimization loop and save data
                    [ obj, ...
                      periods_vs_fills( i_ff_bot, 1 ), ...
                      offsets_vs_fills( i_ff_bot, 1 ), ...
                      directivities_vs_fills( i_ff_bot, 1 ), ...
                      angles_vs_fills( i_ff_bot, 1 ), ...
                      scatter_str_vs_fills( i_ff_bot, 1 ), ...
                      GC_vs_fills{ i_ff_bot, 1 }, ...
                      k_vs_fills( i_ff_bot, 1 ), ...
                      dir_b4_period_vs_fills( i_ff_bot, 1 ) ...
                      ] = ...
                        obj.optimize_period_offset( guess_offset, ...
                                                  fill_top, ...
                                                  fill_bots(i_ff_bot), ...
                                                  guess_period,...
                                                  guessk, ...
                                                  sim_opts, ...
                                                  guess_GC );


                    % update the guess parameters, period, k, offset
                    guessk              = k_vs_fills( i_ff_bot, 1 );
                    guess_period        = periods_vs_fills( i_ff_bot, 1 );
                    guess_GC            = GC_vs_fills{ i_ff_bot, 1 };
                    guess_offset        = offsets_vs_fills( i_ff_bot, 1 );

                else
                    % we're in a waveguide, there's no reason to run
                    % the optimization (and actually the period sweep
                    % bugs out when the fill = 100%)

                    % save dummy 
                    directivities_vs_fills( i_ff_bot, 1 )   = 1;
                    angles_vs_fills( i_ff_bot, 1 )          = 0;
                    scatter_str_vs_fills( i_ff_bot, 1 )     = 0;
                    periods_vs_fills( i_ff_bot, 1 )         = guess_period;
                    offsets_vs_fills( i_ff_bot, 1 )         = 0;
                    k_vs_fills( i_ff_bot, 1 )               = guessk;
                    GC_vs_fills{ i_ff_bot, 1 }              = waveguide;
                    dir_b4_period_vs_fills( i_ff_bot, 1 )   = 1;

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
            n_fill_bots             = length(fill_bots);
            n_fill_top_bot_ratio    = length(fill_top_bot_ratio);

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

                for i_ff_ratio = 2:n_fill_top_bot_ratio
                    % for each top/bottom fill factor ratio

                    fprintf('Fill factor ratio %i of %i, main parfor iteration %i of %i\n', i_ff_ratio, n_fill_top_bot_ratio, i_ff_bot, n_fill_bots);
                        
                    % Optimize period and offset
                    fill_top = fill_top_bot_ratio(i_ff_ratio) * fill_bot;
                    if fill_top < 1
                            % Only run optimization if theres a perturbation 
                            
                        % run optimization loop and save data
                        [ ~, best_period, best_offset, best_directivity, ...
                          best_angle, best_scatter_str, best_GC, best_k, dir_b4_period_vs_fill ] = ...
                            obj.optimize_period_offset( guess_offset, ...
                                                      fill_top, ...
                                                      fill_bot, ...
                                                      guess_period,...
                                                      guessk, ...
                                                      sim_opts, ...
                                                      guess_GC );
                                        
                        periods_vs_fills( i_ff_bot, i_ff_ratio )        = best_period;
                        offsets_vs_fills( i_ff_bot, i_ff_ratio )        = best_offset;
                        directivities_vs_fills( i_ff_bot, i_ff_ratio )  = best_directivity;
                        angles_vs_fills( i_ff_bot, i_ff_ratio )         = best_angle;
                        scatter_str_vs_fills( i_ff_bot, i_ff_ratio )    = best_scatter_str;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }             = best_GC;
                        k_vs_fills( i_ff_bot, i_ff_ratio )              = best_k;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_ratio )  = dir_b4_period_vs_fill;
                
                        % update the guess parameters, period, k, offset
                        guessk              = best_k;
                        guess_period        = best_period;
                        guess_GC            = best_GC;
                        guess_offset        = best_offset;
%                         guessk              = k_vs_fills( i_ff_bot, i_ff_ratio );
%                         guess_period        = periods_vs_fills( i_ff_bot, i_ff_ratio );
%                         guess_GC            = GC_vs_fills{ i_ff_bot, i_ff_ratio };
%                         guess_offset        = offsets_vs_fills( i_ff_bot, i_ff_ratio );
                            
                    else
                        % we're in a waveguide, there's no reason to run
                        % the optimization (and actually the period sweep
                        % bugs out when the fill = 100%)

                        % save dummy 
                        directivities_vs_fills( i_ff_bot, i_ff_ratio )    = 1;
                        angles_vs_fills( i_ff_bot, i_ff_ratio )           = 0;
                        scatter_str_vs_fills( i_ff_bot, i_ff_ratio )      = 0;
                        periods_vs_fills( i_ff_bot, i_ff_ratio )          = guess_period;
                        offsets_vs_fills( i_ff_bot, i_ff_ratio )          = 0;
                        k_vs_fills( i_ff_bot, i_ff_ratio )                = guessk;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }               = waveguide;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_ratio )    = 1;

                    end     % end if fill_top < 1

                end     % end for i_ff_ratio = ...
                
%                 periods_vs_fills( i_ff_bot, 2:end ) = p_v_fill_thisbot;

            end     % end parfor i_ff_bot = ...
            
            % save variables to object
            obj.sweep_variables.directivities_vs_fills  = directivities_vs_fills;
            obj.sweep_variables.angles_vs_fills         = angles_vs_fills;
            obj.sweep_variables.scatter_str_vs_fills    = scatter_str_vs_fills;
            obj.sweep_variables.periods_vs_fills        = periods_vs_fills;
            obj.sweep_variables.offsets_vs_fills        = offsets_vs_fills;
            obj.sweep_variables.k_vs_fills              = k_vs_fills;
%             obj.sweep_variables.GC_vs_fills             = GC_vs_fills;
            obj.sweep_variables.dir_b4_period_vs_fills  = dir_b4_period_vs_fills;
            
            fprintf('Done generating design space\n');
            toc;
            
        end     % end generate_design_space()
        
        
        % -----------------
        % Function generate_design_space_filltopbot()
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
        function obj = generate_design_space_filltopbot( obj, fill_bots, fill_tops, verbose )
            
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
                fill_bots   = fliplr( 0.95:0.025:0.975 );
            end
            if ~exist('fill_tops', 'var')
                % default fill
                fill_tops   = fliplr( 0.95:0.025:0.975 );
            end
            fill_top_bot_ratio  = [];
            guess_offset        = 0;
            
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
            
            % make grating cell, assuming both layers are filled
            waveguide = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.units.name, ...
                                               obj.lambda, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
                                               2*obj.discretization, ...
                                               1.0, 1.0, 0.0 );
            
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
%                 fill_top = fill_top_bot_ratio(1) * fill_bots(i_ff_bot);
                if fill_bots(i_ff_bot) < 1 && fill_tops(1) < 1
                    % Only run optimization if theres a perturbation of
                    % both layers

                    % run optimization loop and save data
                    [ obj, ...
                      periods_vs_fills( i_ff_bot, 1 ), ...
                      offsets_vs_fills( i_ff_bot, 1 ), ...
                      directivities_vs_fills( i_ff_bot, 1 ), ...
                      angles_vs_fills( i_ff_bot, 1 ), ...
                      scatter_str_vs_fills( i_ff_bot, 1 ), ...
                      GC_vs_fills{ i_ff_bot, 1 }, ...
                      k_vs_fills( i_ff_bot, 1 ), ...
                      dir_b4_period_vs_fills( i_ff_bot, 1 ) ...
                      ] = ...
                        obj.optimize_period_offset( guess_offset, ...
                                                  fill_tops(1), ...
                                                  fill_bots(i_ff_bot), ...
                                                  guess_period,...
                                                  guessk, ...
                                                  sim_opts, ...
                                                  guess_GC );


                    % update the guess parameters, period, k, offset
                    guessk              = k_vs_fills( i_ff_bot, 1 );
                    guess_period        = periods_vs_fills( i_ff_bot, 1 );
                    guess_GC            = GC_vs_fills{ i_ff_bot, 1 };
                    guess_offset        = offsets_vs_fills( i_ff_bot, 1 );

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
                            
                        % run optimization loop and save data
                        [ ~, best_period, best_offset, best_directivity, ...
                          best_angle, best_scatter_str, best_GC, best_k, dir_b4_period_vs_fill ] = ...
                            obj.optimize_period_offset( guess_offset, ...
                                                      fill_top, ...
                                                      fill_bot, ...
                                                      guess_period,...
                                                      guessk, ...
                                                      sim_opts, ...
                                                      guess_GC );
                                        
                        periods_vs_fills( i_ff_bot, i_ff_top )        = best_period;
                        offsets_vs_fills( i_ff_bot, i_ff_top )        = best_offset;
                        directivities_vs_fills( i_ff_bot, i_ff_top )  = best_directivity;
                        angles_vs_fills( i_ff_bot, i_ff_top )         = best_angle;
                        scatter_str_vs_fills( i_ff_bot, i_ff_top )    = best_scatter_str;
%                         GC_vs_fills{ i_ff_bot, i_ff_ratio }             = best_GC;
                        k_vs_fills( i_ff_bot, i_ff_top )              = best_k;
                        dir_b4_period_vs_fills( i_ff_bot, i_ff_top )  = dir_b4_period_vs_fill;
                
                        % update the guess parameters, period, k, offset
                        guessk              = best_k;
                        guess_period        = best_period;
                        guess_GC            = best_GC;
                        guess_offset        = best_offset;
%                         guessk              = k_vs_fills( i_ff_bot, i_ff_ratio );
%                         guess_period        = periods_vs_fills( i_ff_bot, i_ff_ratio );
%                         guess_GC            = GC_vs_fills{ i_ff_bot, i_ff_ratio };
%                         guess_offset        = offsets_vs_fills( i_ff_bot, i_ff_ratio );
                            
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

                    end     % end if fill_top < 1

                end     % end for i_ff_top = ...
                
%                 periods_vs_fills( i_ff_bot, 2:end ) = p_v_fill_thisbot;

            end     % end parfor i_ff_bot = ...
            
            % save variables to object
            obj.sweep_variables.directivities_vs_fills  = directivities_vs_fills;
            obj.sweep_variables.angles_vs_fills         = angles_vs_fills;
            obj.sweep_variables.scatter_str_vs_fills    = scatter_str_vs_fills;
            obj.sweep_variables.periods_vs_fills        = periods_vs_fills;
            obj.sweep_variables.offsets_vs_fills        = offsets_vs_fills;
            obj.sweep_variables.k_vs_fills              = k_vs_fills;
%             obj.sweep_variables.GC_vs_fills             = GC_vs_fills;
            obj.sweep_variables.dir_b4_period_vs_fills  = dir_b4_period_vs_fills;
            
            fprintf('Done generating design space\n');
            toc;
            
        end     % end generate_design_space_filltopbot()
        
        
        
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
                    best_GC, best_k, dir_b4_period_vs_fill ] ...
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

                % verbose printing
                if obj.debug_options.verbose == true
                    fprintf('Sweeping offset %i of %i\n', i_offset, length(offset_ratios) );
                end
                
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.units.name, ...
                                               obj.lambda, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
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
            
            % set while loop exit flag
            angle_err_sign_flip = false;

            i_period    = 0;
            while ~angle_err_sign_flip

                i_period = i_period + 1;
                
                % verbose printing
                if obj.debug_options.verbose == true
                    fprintf('Sweeping period %i\n', i_period );
                end
                
                % make grating cell
                GC = obj.h_makeGratingCell( obj.discretization, ...
                                               obj.units.name, ...
                                               obj.lambda, ...
                                               obj.background_index, ...
                                               obj.y_domain_size, ...
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
        
        
        function obj = generate_final_design_gaussian( obj, MFD, input_wg_type, invert_normal_top_bot_ratio_thresh, enforce_min_feat_size_func )
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
            %   invert_normal_top_bot_ratio_thresh
            %       type: double, scalar
            %       desc: OPTIONAL: top/bottom fill factor threshold
            %             currently for my own debugging purposes
            %
            % Sets these fields: (not updated)
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
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % generate x coordinates for the gaussian mode
            % must be large enough to fit mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiber_mode_gaussian(  w0, zvec, xvec,...
                                                    obj.optimal_angle, d0, obj.background_index );
                                              
            % calculate desired scattering strength vs. x
            integral_u      = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
            alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );             % in units 1/m
            alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
            
            % DEBUG plot alpha desired
            figure;
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice();      

            % meshgrid the fills
            [ topbot_ratio_mesh, bot_fills_mesh ] = meshgrid( obj.sweep_variables.fill_top_bot_ratio, obj.sweep_variables.fill_bots );
            
            % set invert/normal threshold
            if ~exist( 'invert_normal_top_bot_ratio_thresh', 'var' )
                invert_normal_top_bot_ratio_thresh = 0.6;
            end

            
            % Select subset of domain to use
            if strcmp( input_wg_type, 'bottom' ) == true
                % input waveguide is the bottom layer
            
                % first narrow down the space to take the datapoints with the
                % maximum directivity per bottom fill factor (so for each row
                % on my design space)

                indx_bottom_space   = obj.sweep_variables.fill_top_bot_ratio < invert_normal_top_bot_ratio_thresh;
                
                % narrow down design space (take top left quadrant)
                % recall fills are sorted in descending order
                topbot_ratios   = topbot_ratio_mesh( :, indx_bottom_space );
                bot_fills       = bot_fills_mesh( :, indx_bottom_space );
                directivities   = obj.sweep_variables.directivities_vs_fills( :, indx_bottom_space );
                angles          = obj.sweep_variables.angles_vs_fills( :, indx_bottom_space );      
                periods         = obj.sweep_variables.periods_vs_fills( :, indx_bottom_space );
                offsets         = obj.sweep_variables.offsets_vs_fills( :, indx_bottom_space );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( :, indx_bottom_space );
                ks              = obj.sweep_variables.k_vs_fills( :, indx_bottom_space );
                
            elseif strcmp( input_wg_type, 'full' ) == true
                % input waveguide is a full waveguide
                
                % first narrow down the space to take the datapoints with the
                % maximum directivity per bottom fill factor (so for each row
                % on my design space)

                indx_full_space   = obj.fill_top_bot_ratio > invert_normal_top_bot_ratio_thresh - 0.01;
                
                % narrow down design space (take top right quadrant)
                topbot_ratios   = topbot_ratio_mesh( :, indx_full_space );
                bot_fills       = bot_fills_mesh( :, indx_full_space );
                directivities   = obj.sweep_variables.directivities_vs_fills( :, indx_full_space );
                angles          = obj.sweep_variables.angles_vs_fills( :, indx_full_space );      
                periods         = obj.sweep_variables.periods_vs_fills( :, indx_full_space );
                offsets         = obj.sweep_variables.offsets_vs_fills( :, indx_full_space );
                scatter_strs    = obj.sweep_variables.scatter_str_vs_fills( :, indx_full_space );
                ks              = obj.sweep_variables.k_vs_fills( :, indx_full_space );
                
            end     % end if strcmp( input_wg_type, 'bottom' )
            
            % for each top fill, pick bottom fill with highest
            % directivity
            [ high_dirs, indx_highest_dir_per_ratio ] = max( directivities, [], 1 );
            % gotta use linear indexing
            indxs                   = sub2ind( size(directivities), indx_highest_dir_per_ratio, 1:size(directivities,2)  );
            angles_high_dir         = angles( indxs );
            periods_high_dir        = periods( indxs );
            offsets_high_dir        = offsets( indxs );
            scatter_strs_high_dir   = scatter_strs( indxs );
            k_high_dir              = ks( indxs );
            topbot_ratio_high_dir   = topbot_ratios( indxs );
            bot_fills_high_dir      = bot_fills( indxs );
%                 GC_high_dir             = GC_vs_fills_inv( indxs );
            

            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;
                
                indices_to_keep = [];
                for ii = 1:length( periods_high_dir )
                    if enforce_min_feat_size_func(  periods_high_dir(ii), ...
                                                    topbot_ratio_high_dir(ii) .* bot_fills_high_dir(ii), ...
                                                    bot_fills_high_dir(ii) ) ...
                                                    == true
                        indices_to_keep(end+1) = ii;
                    end
                end
                
                angles_high_dir         = angles_high_dir(indices_to_keep);
                periods_high_dir        = periods_high_dir(indices_to_keep);
                offsets_high_dir        = offsets_high_dir(indices_to_keep);
                scatter_strs_high_dir   = scatter_strs_high_dir(indices_to_keep);
                k_high_dir              = k_high_dir(indices_to_keep);
                topbot_ratio_high_dir   = topbot_ratio_high_dir(indices_to_keep);
                bot_fills_high_dir      = bot_fills_high_dir(indices_to_keep);
                        
            end

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
            plot( 1:length(high_dirs), 10*log10(high_dirs), '-o' );
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
            
            % match datapoints to desired alpha
            obj = obj.pick_final_datapoints( xvec, alpha_des, high_dirs, ...
                                         bot_fills_high_dir, top_fills_high_dir, ...
                                         offsets_high_dir, periods_high_dir, ...
                                         angles_high_dir, scatter_strs_high_dir, ...
                                         k_high_dir );

%             % now match these data points to the desired alpha
%             % starting point
%             [~, indx_max_alpha] = max( alpha_des );
%             start_alpha_des     = min(scatter_strs_high_dir); %1e-5;
%             [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
%             cur_x               = xvec(indx_x);
%             
%             % final synthesized variables
%             obj.synthesized_design.dir                  = [];
%             obj.synthesized_design.bot_fill             = [];
%             obj.synthesized_design.top_bot_fill_ratio   = [];
%             obj.synthesized_design.top_fill             = [];
%             obj.synthesized_design.period               = [];
%             obj.synthesized_design.offset               = [];
%             obj.synthesized_design.angles               = [];
%             obj.synthesized_design.scatter_str          = [];
%             obj.synthesized_design.k                    = [];
%             obj.synthesized_design.GC                   = {};
%             obj.synthesized_design.des_scatter          = [];
%             
%             
%             % flag for switching to using max scattering strength
%             saturate_scatter_str_to_max = false;
%  
%             ii = 1;
%             while cur_x < xvec(end)
%                 % build grating one cell at a time
%                 
%                 % pick design with scattering strength closest to desired
%                 % alpha
%                 des_scatter                 = alpha_des(indx_x);                                        % desired alpha
%                 if des_scatter  > max( scatter_strs_high_dir )
%                     % desired scattering strength too high, gotta saturate
%                     saturate_scatter_str_to_max = true;
%                 end
%                 if ~saturate_scatter_str_to_max
%                     [~, indx_closest_scatter]   = min( abs(scatter_strs_high_dir - des_scatter) );          % index of closest scatter design 
%                 else
%                     [~, indx_closest_scatter]   = max( scatter_strs_high_dir );                             % saturate to max
%                 end
%                 
%                 % save parameters
%                 obj.synthesized_design.dir(ii)                   = high_dirs( indx_closest_scatter );
%                 obj.synthesized_design.bot_fill(ii)              = bot_fills_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.top_bot_fill_ratio(ii)    = topbot_ratio_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.top_fill(ii)              = topbot_ratio_high_dir( indx_closest_scatter ) * bot_fills_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.offset(ii)                = offsets_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.period(ii)                = periods_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.angles(ii)                = angles_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.scatter_str(ii)           = scatter_strs_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.k(ii)                     = k_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.des_scatter(ii)           = des_scatter;
%                 
%                 obj.synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
%                                                        obj.discretization, ...
%                                                        obj.units.name, ...
%                                                        obj.lambda, ...
%                                                        obj.background_index, ...
%                                                        obj.y_domain_size, ...
%                                                        obj.synthesized_design.period(ii), ...
%                                                        obj.synthesized_design.top_fill(ii), ...
%                                                        obj.synthesized_design.bot_fill(ii), ...
%                                                        obj.synthesized_design.offset(ii)/obj.synthesized_design.period(ii) );
% 
%                 
%                 % move onto next
%                 cur_x       = cur_x + obj.synthesized_design.period(ii);
%                 [~, indx_x] = min( abs(xvec - cur_x) );
%                 cur_x       = xvec( indx_x );
%                 ii          = ii + 1;
%                 
%             end     % end for ii = 1:ncells
            
            
            % build final index distribution
            obj = obj.build_final_index();
%             obj.synthesized_design.N = [];
%             for ii = 1:length(obj.synthesized_design.GC)
%                
%                 GC                          = obj.synthesized_design.GC{ii};
%                 obj.synthesized_design.N    = [ obj.synthesized_design.N, GC.N ];
%                 
%             end
            
            % coordinates of index distribution
            obj.synthesized_design.x_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,2)-1 ) );
            obj.synthesized_design.y_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,1)-1 ) );
            
            
        end     % end function generate_final_design_gaussian()
        
        
        function obj = generate_final_design_gaussian_topbot( obj, MFD, input_wg_type, enforce_min_feat_size_func )
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
            %   enforce_min_feat_size_func
            %       type: function handle
            %       desc: OPTIONAL INPUT
            %             A function that user makes which enforces min.
            %             feat size
            %             See the example function at the bottom of this
            %             file
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
            

            % save input waveguide type
            obj.synthesized_design.input_wg_type = input_wg_type;
            
            % generate x coordinates for the gaussian mode
            % must be large enough to fit mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiber_mode_gaussian(  w0, zvec, xvec,...
                                                    obj.optimal_angle, d0, obj.background_index );
                                              
            % calculate desired scattering strength vs. x
            integral_u      = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
            alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );             % in units 1/m
            alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
            
            % DEBUG plot alpha desired
            figure;
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice();      

            % meshgrid the fills (dimensions are bot vs. top, I believe)
            [ top_fills_mesh, bot_fills_mesh ] = meshgrid( obj.sweep_variables.fill_tops, obj.sweep_variables.fill_bots );

            % find where the max scattering strength is
            [~, max_indx] = max( obj.sweep_variables.scatter_str_vs_fills(:) );
            [ i_max_bot_scat, i_max_top_scat ] = ind2sub( size(obj.sweep_variables.scatter_str_vs_fills), max_indx );
            
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
                
            end     % end if strcmp( input_wg_type, 'bottom' )
            
            % for each top fill, pick bottom fill with highest
            % directivity
            [ high_dirs, indx_highest_dir_per_topfill ] = max( directivities, [], 1 );
            % gotta use linear indexing
            indxs                   = sub2ind( size(directivities), indx_highest_dir_per_topfill, 1:size(directivities,2)  );
            angles_high_dir         = angles( indxs );
            periods_high_dir        = periods( indxs );
            offsets_high_dir        = offsets( indxs );
            scatter_strs_high_dir   = scatter_strs( indxs );
            k_high_dir              = ks( indxs );
            top_fills_high_dir      = top_fills( indxs );
            bot_fills_high_dir      = bot_fills( indxs );
%                 GC_high_dir             = GC_vs_fills_inv( indxs );

            % enforce min feature size
            obj.synthesized_design.use_min_feat_size = false;                   % default to false
            if exist( 'enforce_min_feat_size_func', 'var' )
                % loop through each cell and discard any that violate
                % feature size rules
                obj.synthesized_design.use_min_feat_size = true;
                
                indices_to_keep = [];
                for ii = 1:length( periods_high_dir )
                    if enforce_min_feat_size_func( periods_high_dir(ii), top_fills_high_dir(ii), bot_fills_high_dir(ii) ) == true
                        indices_to_keep(end+1) = ii;
                    end
                end
                
                angles_high_dir         = angles_high_dir(indices_to_keep);
                periods_high_dir        = periods_high_dir(indices_to_keep);
                offsets_high_dir        = offsets_high_dir(indices_to_keep);
                scatter_strs_high_dir   = scatter_strs_high_dir(indices_to_keep);
                k_high_dir              = k_high_dir(indices_to_keep);
                top_fills_high_dir      = top_fills_high_dir(indices_to_keep);
                bot_fills_high_dir      = bot_fills_high_dir(indices_to_keep);
                        
            end
                
            
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
            plot( 1:length(top_fills_high_dir), top_fills_high_dir, '-o' );
            title('chosen datapoints, top fill');
            makeFigureNice();
            % bottom fill
            figure;
            plot( 1:length(bot_fills_high_dir), bot_fills_high_dir, '-o' );
            title('chosen datapoints, bottom fill');
            makeFigureNice();
            % directivity
            figure;
            plot( 1:length(high_dirs), 10*log10(high_dirs), '-o' );
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
            
            % match data points to desired alpha
            obj = obj.pick_final_datapoints( xvec, alpha_des, high_dirs, ...
                                         bot_fills_high_dir, top_fills_high_dir, ...
                                         offsets_high_dir, periods_high_dir, ...
                                         angles_high_dir, scatter_strs_high_dir, ...
                                         k_high_dir );

%             % now match these data points to the desired alpha
%             % starting point
%             [~, indx_max_alpha] = max( alpha_des );
%             start_alpha_des     = min(scatter_strs_high_dir); %1e-5;
%             [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
%             cur_x               = xvec(indx_x);
%             
%             % final synthesized variables
%             obj.synthesized_design.dir                  = [];
%             obj.synthesized_design.bot_fill             = [];
%             obj.synthesized_design.top_bot_fill_ratio   = [];
%             obj.synthesized_design.top_fill             = [];
%             obj.synthesized_design.period               = [];
%             obj.synthesized_design.offset               = [];
%             obj.synthesized_design.angles               = [];
%             obj.synthesized_design.scatter_str          = [];
%             obj.synthesized_design.k                    = [];
%             obj.synthesized_design.GC                   = {};
%             obj.synthesized_design.des_scatter          = [];
%             
%             
%             % flag for switching to using max scattering strength
%             saturate_scatter_str_to_max = false;
%  
%             ii = 1;
%             while cur_x < xvec(end)
%                 % build grating one cell at a time
%                 
%                 % pick design with scattering strength closest to desired
%                 % alpha
%                 des_scatter = alpha_des(indx_x);                            % desired alpha
%                 if des_scatter  > max( scatter_strs_high_dir )
%                     % desired scattering strength too high, gotta saturate
%                     saturate_scatter_str_to_max = true;
%                 end
%                 if ~saturate_scatter_str_to_max
%                     [~, indx_closest_scatter]   = min( abs(scatter_strs_high_dir - des_scatter) );          % index of closest scatter design 
%                 else
%                     [~, indx_closest_scatter]   = max( scatter_strs_high_dir );                             % saturate to max
%                 end
%                 
%                 % save parameters
%                 obj.synthesized_design.dir(ii)                   = high_dirs( indx_closest_scatter );
%                 obj.synthesized_design.bot_fill(ii)              = bot_fills_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.top_fill(ii)              = top_fills_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.top_bot_fill_ratio(ii)    = top_fills_high_dir( indx_closest_scatter ) ./ bot_fills_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.offset(ii)                = offsets_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.period(ii)                = periods_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.angles(ii)                = angles_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.scatter_str(ii)           = scatter_strs_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.k(ii)                     = k_high_dir( indx_closest_scatter );
%                 obj.synthesized_design.des_scatter(ii)           = des_scatter;
%                 
%                 obj.synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
%                                                        obj.discretization, ...
%                                                        obj.units.name, ...
%                                                        obj.lambda, ...
%                                                        obj.background_index, ...
%                                                        obj.y_domain_size, ...
%                                                        obj.synthesized_design.period(ii), ...
%                                                        obj.synthesized_design.top_fill(ii), ...
%                                                        obj.synthesized_design.bot_fill(ii), ...
%                                                        obj.synthesized_design.offset(ii)/obj.synthesized_design.period(ii) );
% 
%                 
%                 % move onto next
%                 cur_x       = cur_x + obj.synthesized_design.period(ii);
%                 [~, indx_x] = min( abs(xvec - cur_x) );
%                 cur_x       = xvec( indx_x );
%                 ii          = ii + 1;
%                 
%             end     % end for ii = 1:ncells
            
            
            % build final index distribution
            obj = obj.build_final_index();
%             obj.synthesized_design.N = [];
%             for ii = 1:length(obj.synthesized_design.GC)
%                
%                 GC                          = obj.synthesized_design.GC{ii};
%                 obj.synthesized_design.N    = [ obj.synthesized_design.N, GC.N ];
%                 
%             end
            
            % coordinates of index distribution
            obj.synthesized_design.x_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,2)-1 ) );
            obj.synthesized_design.y_coords = obj.discretization*( 0:1:( size(obj.synthesized_design.N,1)-1 ) );
            
            
        end     % end function generate_final_design_gaussian_topbot()
        
        
        function obj = generate_final_design_uniform( obj, MFD )
            % purpose of this function right now is just for uniform
            % gratings for CLO
            
            % save input waveguide type
            obj.synthesized_design.input_wg_type = 'bottom';
            
            
            
        end
        
        
        function [ obj, alpha_des ] = calculate_desired_scattering( obj )
            % Calculates desired scattering profile for a Gaussian field
            % with the given MFD
           
            % generate x coordinates for the gaussian mode
            % must be large enough to fit mode
            xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
            xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
            
            % generate a fiber gaussian mode
            w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
            zvec        = 0;                                                            % this is unused
            d0          = 0;                                                            % take slice at waist
            [obj, u]    = obj.fiber_mode_gaussian(  w0, zvec, xvec,...
                                                    obj.optimal_angle, d0, obj.background_index );
                                              
            % calculate desired scattering strength vs. x
            integral_u      = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
            alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );             % in units 1/m
            alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
            
            % DEBUG plot alpha desired
            figure;
            plot( xvec, alpha_des );
            xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
            title('DEBUG scattering strength for gaussian');
            makeFigureNice(); 
            
        end
        
        
        function obj = pick_final_datapoints( obj, xvec, alpha_des, high_dirs, ...
                                              bot_fills_high_dir, top_fills_high_dir, ...
                                              offsets_high_dir, periods_high_dir, ...
                                              angles_high_dir, scatter_strs_high_dir, ...
                                              k_high_dir )
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
            
            % now match these data points to the desired alpha
            % starting point
            [~, indx_max_alpha] = max( alpha_des );
            start_alpha_des     = min(scatter_strs_high_dir); %1e-5;
            [~, indx_x]         = min(abs( alpha_des(1:indx_max_alpha) - start_alpha_des ) );
            cur_x               = xvec(indx_x);
            
            % final synthesized variables
            obj.synthesized_design.dir                  = [];
            obj.synthesized_design.bot_fill             = [];
            obj.synthesized_design.top_bot_fill_ratio   = [];
            obj.synthesized_design.top_fill             = [];
            obj.synthesized_design.period               = [];
            obj.synthesized_design.offset               = [];
            obj.synthesized_design.angles               = [];
            obj.synthesized_design.scatter_str          = [];
            obj.synthesized_design.k                    = [];
            obj.synthesized_design.GC                   = {};
            obj.synthesized_design.des_scatter          = [];
            
            
            % flag for switching to using max scattering strength
            saturate_scatter_str_to_max = false;
 
            ii = 1;
            while cur_x < xvec(end)
                % build grating one cell at a time
                
                % pick design with scattering strength closest to desired
                % alpha
                des_scatter = alpha_des(indx_x);                            % desired alpha
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
                obj.synthesized_design.dir(ii)                   = high_dirs( indx_closest_scatter );
                obj.synthesized_design.bot_fill(ii)              = bot_fills_high_dir( indx_closest_scatter );
                obj.synthesized_design.top_fill(ii)              = top_fills_high_dir( indx_closest_scatter );
                obj.synthesized_design.top_bot_fill_ratio(ii)    = top_fills_high_dir( indx_closest_scatter ) ./ bot_fills_high_dir( indx_closest_scatter );
                obj.synthesized_design.offset(ii)                = offsets_high_dir( indx_closest_scatter );
                obj.synthesized_design.period(ii)                = periods_high_dir( indx_closest_scatter );
                obj.synthesized_design.angles(ii)                = angles_high_dir( indx_closest_scatter );
                obj.synthesized_design.scatter_str(ii)           = scatter_strs_high_dir( indx_closest_scatter );
                obj.synthesized_design.k(ii)                     = k_high_dir( indx_closest_scatter );
                obj.synthesized_design.des_scatter(ii)           = des_scatter;
                
                obj.synthesized_design.GC{ii} = obj.h_makeGratingCell(    ...
                                                       obj.discretization, ...
                                                       obj.units.name, ...
                                                       obj.lambda, ...
                                                       obj.background_index, ...
                                                       obj.y_domain_size, ...
                                                       obj.synthesized_design.period(ii), ...
                                                       obj.synthesized_design.top_fill(ii), ...
                                                       obj.synthesized_design.bot_fill(ii), ...
                                                       obj.synthesized_design.offset(ii)/obj.synthesized_design.period(ii) );

                
                % move onto next
                cur_x       = cur_x + obj.synthesized_design.period(ii);
                [~, indx_x] = min( abs(xvec - cur_x) );
                cur_x       = xvec( indx_x );
                ii          = ii + 1;
                
            end     % end for ii = 1:ncells
            
        end     % end pick_final_datapoints()
        
        
        function obj = build_final_index( obj )
        % build final index distribution
        % to be run after a design has been synthesized
        
            obj.synthesized_design.N = [];
            for ii = 1:length(obj.synthesized_design.GC)

                GC                          = obj.synthesized_design.GC{ii};
                obj.synthesized_design.N    = [ obj.synthesized_design.N, GC.N ];

            end
            
        end     % end build_final_index()
        
%         function obj = generateFinalDesignGaussian_old(obj, MFD, input_wg_type)
%             % OLD VERSION
%             % function for generating the final synthesized design
%             % parameters
%             %
%             % Inputs:
%             %   MFD
%             %       type: double, scalar
%             %       desc: mode field diameter, in units 'units'
%             %   input_wg_type
%             %       type: string
%             %       desc: 'bottom' for cSi only or 'full' for both cSi and
%             %             pSi
%             %
%             % Sets these fields:
%             %   obj.dir_synth                   = [];
%             %   obj.bot_fill_synth              = [];
%             %   obj.top_bot_fill_ratio_synth    = [];
%             %   obj.period_synth                = [];
%             %   obj.offset_synth                = [];
%             %   obj.angles_synth                = [];
%             %   obj.scatter_str_synth           = [];
%             %   obj.k_synth                     = [];
%             %   obj.GC_synth                    = {};
%             %   obj.des_scatter_norm            = [];
%             
% 
%             % save input waveguide type
%             obj.input_wg_type = input_wg_type;
%             
%             % generate x coordinates for the gaussian mode
%             % must be large enough to fit mode
%             xvec            = 0 : obj.discretization : MFD*4 - obj.discretization;
%             xvec            = xvec - xvec(round(end/2));                                % shift origin over to middle
%             
%             % generate a fiber gaussian mode
%             w0          = MFD/2;                                                        % not sure if this is the proper exact relationship
%             zvec        = 0;                                                            % this is unused
%             d0          = 0;                                                            % take slice at waist
%             [obj, u]    = obj.fiberModeGaussian(    w0, zvec, xvec,...
%                                                     obj.optimal_angle, d0, obj.background_index );
%             
%             % calculate desired scattering strength vs. x
%             integral_u      = cumsum( abs(u).^2 ) * obj.discretization * obj.units.scale;
%             alpha_des       = (1/2)*( abs(u).^2 ) ./ ( 1 + 1e-9 - integral_u );  % in units 1/m
% %             alpha_des(end)  = 0;                                                            % for stability
%             alpha_des       = alpha_des * obj.units.scale;                                  % in units 1/units
%             
%             
%             % DEBUG plot alpha desired
%             figure;
%             plot( xvec, alpha_des );
%             xlabel(['x (' obj.units.name ')']); ylabel( ['\alpha (1/' obj.units.name ')'] );
%             title('DEBUG scattering strength for gaussian');
%             makeFigureNice();
% %             
% 
%             % meshgrid the fills
%             [ topbot_ratio_mesh, bot_fills_mesh ] = meshgrid( obj.fill_top_bot_ratio, obj.fill_bots );
%             
%             % set invert/normal threshold
%             invert_normal_top_bot_ratio_thresh = 0.6;
% 
%             if strcmp( input_wg_type, 'bottom' ) == true
%                 % Inverted design
%             
%                 % first narrow down the space to take the datapoints with the
%                 % maximum directivity per bottom fill factor (so for each row
%                 % on my design space)
% 
%                 indx_bottom_space   = obj.fill_top_bot_ratio < invert_normal_top_bot_ratio_thresh;
% 
%                 % grab variables (remember dimensions are bot fill x top bot
%                 % ratio)
%                 topbot_ratio_vs_fills_inv   = topbot_ratio_mesh( :, indx_bottom_space );
%                 bot_fills_vs_fills_inv      = bot_fills_mesh( :, indx_bottom_space );
%                 directivities_vs_fills_inv  = obj.directivities_vs_fills( :, indx_bottom_space );
%                 angles_vs_fills_inv         = obj.angles_vs_fills( :, indx_bottom_space );      
%                 periods_vs_fills_inv        = obj.periods_vs_fills( :, indx_bottom_space );
%                 offsets_vs_fills_inv        = obj.offsets_vs_fills( :, indx_bottom_space );
%                 scatter_str_vs_fills_inv    = obj.scatter_str_vs_fills( :, indx_bottom_space );
%                 k_vs_fills_inv              = obj.k_vs_fills( :, indx_bottom_space );
%                 GC_vs_fills_inv             = obj.GC_vs_fills( :, indx_bottom_space );
%                 
%                 % remove datapoints where the angle deviates beyond some
%                 % angle, by artifically setting the directivity to be very
%                 % low
%                 directivities_vs_fills_inv( abs( angles_vs_fills_inv - obj.optimal_angle ) > 1 ) = 1e-6;
% 
% 
%                 % for each top/bottom ratio, pick bottom fill with highest
%                 % directivity
%                 [ highest_dir_per_ratio, indx_highest_dir_per_ratio ] = max( directivities_vs_fills_inv, [], 1 );
%                 % gotta use linear indexing
%                 indxs                   = sub2ind( size(angles_vs_fills_inv), indx_highest_dir_per_ratio, 1:size(angles_vs_fills_inv,2)  );
%                 angles_high_dir         = angles_vs_fills_inv( indxs );
%                 periods_high_dir        = periods_vs_fills_inv( indxs );
%                 offsets_high_dir        = offsets_vs_fills_inv( indxs );
%                 scatter_strs_high_dir   = scatter_str_vs_fills_inv( indxs );
%                 k_high_dir              = k_vs_fills_inv( indxs );
%                 topbot_ratio_high_dir   = topbot_ratio_vs_fills_inv( indxs );
%                 bot_fills_high_dir      = bot_fills_vs_fills_inv( indxs );
%                 GC_high_dir             = GC_vs_fills_inv( indxs );
%                 dir_high                = highest_dir_per_ratio;
%                 
%             elseif strcmp( input_wg_type, 'full' ) == true
%                 % normal design
%                 
%                 % first narrow down the space to take the datapoints with the
%                 % maximum directivity per bottom fill factor (so for each row
%                 % on my design space)
% 
%                 indx_full_space   = obj.fill_top_bot_ratio > invert_normal_top_bot_ratio_thresh - 0.01;
% 
%                 % grab variables (remember dimensions are bot fill x top bot
%                 % ratio)
%                 topbot_ratio_vs_fills_full   = topbot_ratio_mesh( :, indx_full_space );
%                 bot_fills_vs_fills_full      = bot_fills_mesh( :, indx_full_space );
%                 directivities_vs_fills_full  = obj.directivities_vs_fills( :, indx_full_space );
%                 angles_vs_fills_full         = obj.angles_vs_fills( :, indx_full_space );      
%                 periods_vs_fills_full        = obj.periods_vs_fills( :, indx_full_space );
%                 offsets_vs_fills_full        = obj.offsets_vs_fills( :, indx_full_space );
%                 scatter_str_vs_fills_full    = obj.scatter_str_vs_fills( :, indx_full_space );
%                 k_vs_fills_full              = obj.k_vs_fills( :, indx_full_space );
%                 GC_vs_fills_full             = obj.GC_vs_fills( :, indx_full_space );
%                 
%                 % remove datapoints where the angle deviates beyond some
%                 % angle, by artifically setting the directivity to be very
%                 % low
%                 directivities_vs_fills_full( abs( angles_vs_fills_full - obj.optimal_angle ) > 0.5 ) = 1e-6;
%                 
%                 % for each bottom fill, pick the datapoint with the highest
%                 % directivity
%                 [ highest_dir_per_bot, indx_highest_dir_per_bot ] = max( directivities_vs_fills_full, [], 2 );
%                 % gotta use linear indexing
%                 indxs                   = sub2ind( size(angles_vs_fills_full), 1:size(angles_vs_fills_full,1), indx_highest_dir_per_bot.' );
%                 angles_high_dir         = angles_vs_fills_full( indxs );
%                 periods_high_dir        = periods_vs_fills_full( indxs );
%                 offsets_high_dir        = offsets_vs_fills_full( indxs );
%                 scatter_strs_high_dir   = scatter_str_vs_fills_full( indxs );
%                 k_high_dir              = k_vs_fills_full( indxs );
%                 topbot_ratio_high_dir   = topbot_ratio_vs_fills_full( indxs );
%                 bot_fills_high_dir      = bot_fills_vs_fills_full( indxs );
%                 GC_high_dir             = GC_vs_fills_full( indxs );
%                 dir_high                = highest_dir_per_bot;
%                 
%             end     % end if strcmp( input_wg_type, 'bottom' )
%             
%             % DEBUG plot the resulting picked out datapoints
%             % angles
%             figure;
%             plot( 1:length(angles_high_dir), angles_high_dir, '-o' );
%             title('chosen datapoints, angles');
%             makeFigureNice();
%             % periods
%             figure;
%             plot( 1:length(periods_high_dir), periods_high_dir, '-o' );
%             title('chosen datapoints, periods');
%             makeFigureNice();
%             % offsets
%             figure;
%             plot( 1:length(offsets_high_dir), offsets_high_dir, '-o' );
%             title('chosen datapoints, offsets');
%             makeFigureNice();
%             % scattering strengths
%             figure;
%             plot( 1:length(scatter_strs_high_dir), scatter_strs_high_dir, '-o' );
%             title('chosen datapoints, scattering strengths');
%             makeFigureNice();
%             % k real
%             figure;
%             plot( 1:length(k_high_dir), real(k_high_dir), '-o' );
%             title('chosen datapoints, k real');
%             makeFigureNice();
%             % k imag
%             figure;
%             plot( 1:length(k_high_dir), imag(k_high_dir), '-o' );
%             title('chosen datapoints, k imaginary');
%             makeFigureNice();
%             % top bottom ratio
%             figure;
%             plot( 1:length(topbot_ratio_high_dir), topbot_ratio_high_dir, '-o' );
%             title('chosen datapoints, top/bottom ratio');
%             makeFigureNice();
%             % bottom fill
%             figure;
%             plot( 1:length(bot_fills_high_dir), bot_fills_high_dir, '-o' );
%             title('chosen datapoints, bottom fill');
%             makeFigureNice();
%             % directivity
%             figure;
%             plot( 1:length(dir_high), 10*log10(dir_high), '-o' );
%             title('chosen datapoints, directivity (dB)');
%             makeFigureNice();
%             
% %             % DEBUG plot bot fills vs topbot ratio
% %             figure;
% %             plot( topbot_ratio_high_dir, bot_fills_high_dir, '-o' );
% %             xlabel('top/bottom ratio'); ylabel('bottom fill');
% %             title('DEBUG bottom fill vs top/bottom ratio');
% %             makeFigureNice();
% %             
% %             % DEBUG plot bot fills vs top fills
% %             figure;
% %             plot( bot_fills_high_dir, topbot_ratio_high_dir .* bot_fills_high_dir, '-o' );
% %             xlabel('bottom fill'); ylabel('top fill');
% %             title('DEBUG bottom fill vs top/bottom ratio');
% %             makeFigureNice();
% %             
% %             p = polyfit( topbot_ratio_high_dir, bot_fills_high_dir, 2 )
% %             
% %             x = 0:0.01:1;
% %             y = -(x.^3)/2 + x;
% %             figure;
% %             plot(x, y);
% %             xlabel('bot'); ylabel('top');
%             
%             % now match these data points to the desired alpha
%             % starting point
%             start_alpha_des     = 1e-5;
%             [~, indx_x]         = min(abs( alpha_des - start_alpha_des ) );
%             cur_x               = xvec(indx_x);
%             
%             % final synthesized variables
%             obj.synthesized_design.dir_synth                   = [];
%             obj.synthesized_design.bot_fill_synth              = [];
%             obj.synthesized_design.top_bot_fill_ratio_synth    = [];
%             obj.synthesized_design.period_synth                = [];
%             obj.synthesized_design.offset_synth                = [];
%             obj.synthesized_design.angles_synth                = [];
%             obj.synthesized_design.scatter_str_synth           = [];
%             obj.synthesized_design.k_synth                     = [];
%             obj.synthesized_design.GC_synth                    = {};
%             obj.synthesized_design.des_scatter_synth           = [];
%             
%             % flag for switching to using max scattering strength
%             saturate_scatter_str_to_max = false;
%  
%             ii = 1;
%             while cur_x < xvec(end)
%                 % build grating one cell at a time
%                 
%                 % pick design with scattering strength closest to desired
%                 % alpha
%                 des_scatter                 = alpha_des(indx_x);                                        % desired alpha
%                 if des_scatter  > max( scatter_strs_high_dir )
%                     % desired scattering strength too high, gotta saturate
%                     saturate_scatter_str_to_max = true;
%                 end
%                 if ~saturate_scatter_str_to_max
%                     [~, indx_closest_scatter]   = min( abs(scatter_strs_high_dir - des_scatter) );          % index of closest scatter design 
%                 else
%                     [~, indx_closest_scatter]   = max( scatter_strs_high_dir );                             % saturate to max
%                 end
%                 
%                 % save parameters
%                 obj.dir_synth(ii)                   = dir_high( indx_closest_scatter );
%                 obj.bot_fill_synth(ii)              = bot_fills_high_dir( indx_closest_scatter );
%                 obj.top_bot_fill_ratio_synth(ii)    = topbot_ratio_high_dir( indx_closest_scatter );
%                 obj.offset_synth(ii)                = offsets_high_dir( indx_closest_scatter );
%                 obj.period_synth(ii)                = periods_high_dir( indx_closest_scatter );
%                 obj.angles_synth(ii)                = angles_high_dir( indx_closest_scatter );
%                 obj.scatter_str_synth(ii)           = scatter_strs_high_dir( indx_closest_scatter );
%                 obj.k_synth(ii)                     = k_high_dir( indx_closest_scatter );
% %                 obj.GC_synth{ii}                    = GC_high_dir{ indx_closest_scatter };
%                 obj.des_scatter_synth(ii)           = des_scatter;
%                 
%                 % move onto next
%                 cur_x       = cur_x + obj.period_synth(ii);
%                 [~, indx_x] = min( abs(xvec - cur_x) );
%                 cur_x       = xvec( indx_x );
%                 ii          = ii + 1;
%                 
%             end     % end for ii = 1:ncells
%             
%             
%         end     % end function generateFinalDesignGaussian_old()
        
        
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
            figure;
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
















































