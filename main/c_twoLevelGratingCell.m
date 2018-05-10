classdef c_twoLevelGratingCell
% authors: bohan zhang
% 
% encapsulates a 2 level grating cell simulation
% adapted from mark's code
%
% Notes:
%   - propagation direction coordinate is 'x'. 
%   - transverse (in plane) direction coordinate is 'y'.
%   - transverse (out of plane) direction is 'z'.
%
% Inputs to constructor:
%   Inputs are name-value pairs:
%   'discretization'
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y, in units of 'units'
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
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
%   'num_cells'
%       type: integer, scalar
%       desc: OPTIONAL: pick the number of grating cells to repeat when
%             drawing/calculating Ez
%
% 
% Example usage:
%   GC = c_twoLevelGratingCell( 'discretization', discretization, ...
%                               'units', 'nm', ...
%                               'lambda', 1550, ...
%                               'domain_size', [ 2000, 800], ...
%                               'background_index', 1.0, ...
%                               'numcells', 10 );
    
    properties

        N;       % index profile
        dx;      % discretization in x (dir of propagation)
        dy;      % discretization in y (transverse direction)

        % -----------
        % new properties
        units;          % name and scaling of spatial units, supports 'm', 'mm', 'um', 'nm'
        lambda;
        domain_size;    % [ max_y, max_x ], where y = transverse and x = direction of propagation
        x_coords;       % vector of x coordinates (dir of prop)
        y_coords;       % vector of y coords (transverse dir)
        wg_min_y;       % bottom position of wg
        wg_max_y;       % top position of wg
        
        
        % struct that holds simulation options
        % current options are:
        %   'num_modes'     - number of modes
        %   'BC'            - boundary conditions, 0 for pec 1 for pmc i think
        %   'pml_options'   - pml options, see complexk solver for details
        %                     PML_options(1): PML in y direction (yes=1 or no=0)
        %                     PML_options(2): length of PML layer in nm
        %                     PML_options(3): strength of PML in the complex plane
        %                     PML_options(4): PML polynomial order (1, 2, 3...)
        sim_opts;
        
        % mode characteristics
        k;                  % complex k, units rad/'units'
        Phi;                % field envelope
        numcells;           % number of cells repeated in E_z
        E_z                 % saves n_periods of the field, Phi(x,y)*exp(jkx) (z stands for z polarization)
        k_vs_mode;          % k vs mode #
        Phi_vs_mode;        % Phi vs mode #, dimensions y vs. x vs mode #
        chosen_mode_num;    % which mode was chosen
        E_z_for_overlap;    % Field for mode overlapping. One unit cell with only real(k) in the phase
        
        % radiation parameters
        directivity;        % ratio of power radiated up/power radiated down
        max_angle_up;       % Angle of maximum upwards output radiation
        max_angle_down;     % Angle of maximum downwards output radiation
        P_rad_down;         % saves power radiated down
        P_rad_up;           % saves power radiated up
        P_in;               % saves power entering the wg
        alpha_up;           % fractional radiative loss/unit length upwards
        alpha_down;         % fractional radiative loss/unit length downwards
        Sx;                 % X component of poynting vector (dims y vs. x(2:end-1))
        P_per_x_slice;      % x traveling power per x slice vs. x(2:end-1)
        Sy;                 % Y component of poynting vector (dims y(2:end-1) vs x)
        P_per_y_slice;      % Y traveling (upwards) power per y slice vs. y(2:end-1)
        
        % DEBUG struct for holding temporary values that are useful during
        % debugging
        %   Current fields: k_all, phi_all, unguided_power, guided_power,
        %                   p_rad_up, Sx_up, Sx_down, Sy_up, Sy_down, Sx_in
        %                   P_rad_up_onecell, P_rad_down_onecell,
        %                   Sx, Sy, P_per_y_slice, P_per_x_slice
        %   Some of the above fields may have been removed.
        debug;
        
    end
    
    methods
        
        function obj = c_twoLevelGratingCell(varargin)
            % Constructor
            % See top of class definition for input documentation

%             % add paths
%             % TEMPORARY -> path to complex k mode solver
%             addpath('..');
            
            % parse inputs
            inputs = {  'discretization',   'none', ...
                        'units',            'um',   ...
                        'lambda',           'none', ...
                        'background_index', 1.0,    ...
                        'domain_size',      'none', ...
                        'numcells',         10 ...
                     }; 
                 
            p = f_parse_varargin( inputs, varargin{:} );
            
            % sort the inputs
            if length(p.discretization) == 1
                obj.dx          = p.discretization;
                obj.dy          = p.discretization;
            elseif length(p.discretization) == 2
                obj.dx          = p.discretization(2);
                obj.dy          = p.discretization(1);
                if abs(obj.dx - obj.dy) >= 1e-5
                    warning('dx and dy are different. Bloch solver does not currently support different discretizations for different axes');
                end
            else
                error('Input ''discretization'' must either be a scaler or a 1x2 vector');
            end
            
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
            
            obj.lambda = p.lambda;
            
            
            % new version of creating backgrouund dielectric that rounds
            % the dimensions
            obj.domain_size = p.domain_size;
            nx              = round( obj.domain_size(2)/obj.dx );                                           % number of x samples
            obj.x_coords    = 0 : obj.dx : ( (nx-1) * obj.dx );
            ny              = round( obj.domain_size(1)/obj.dy );                                           % number of y samples
            obj.y_coords    = 0 : obj.dy : ( (ny-1) * obj.dy );
            obj.N           = p.background_index * ones( length( obj.y_coords ), length( obj.x_coords ) );  % dimensions of y vs. x
            
            
            % check discretization fits in integer amt
            if abs(obj.x_coords(end) - (obj.domain_size(2) - obj.dx)) >= 1e-10
                % the discretization doesn't fit into x
                fprintf('Warning: dx doesn''t fit integer times into the x domain size.\n');
                
                % override the domain size
                fprintf('Overriding domain. Old domain x size was: %f. New domain x size is: %f\n\n',  obj.domain_size(2), obj.x_coords(end) + obj.dx);
                obj.domain_size(2) = obj.x_coords(end) + obj.dx;
                obj.domain_size(2) = obj.dx * round( obj.domain_size(2)/obj.dx );   % round
            end
            if abs(obj.y_coords(end) - (obj.domain_size(1) - obj.dy)) >= 1e-10
                % the discretization doesn't fit into x
                fprintf('Warning: dy doesn''t fit integer times into the y domain size.\n');
                
                % override the domain size
                fprintf('Overriding domain. Old domain y size was: %f. New domain y size is: %f\n\n', obj.domain_size(1), obj.y_coords(end) + obj.dy);
                obj.domain_size(1) = obj.y_coords(end) + obj.dy;
                obj.domain_size(1) = obj.dy * round( obj.domain_size(1)/obj.dy );   % round
            end
            
            % set number of cells
            obj.numcells = p.numcells;
            
        end     % end constructor
        
        
        function obj = addLayer( obj, min_y, height_y, index )
            % Draws a horizontal layer of dielectric
            %
            % Inputs:
            %   min_y
            %       Desc: scalar double, minimum y
            %   height_y
            %       Desc: scalar double, height/thickness of layer
            %   index
            %       Desc: scalar double, index of refraction of layer
            
            y = obj.y_coords;
            
            % fill in the layer
            obj.N( y > (min_y - obj.dy/2) & y <= (min_y + height_y - obj.dy/2), : ) = index; 
            
        end     % end function addLayer()
        
        
        function obj = addRect( obj, min_x, min_y, width_x, height_y, index )
            % Draws a rectangle
            %
            % Inputs:
            %   min_x
            %       Desc: scalar double, left edge of rectangle
            %   min_y
            %       Desc: scalar double, bottom edge of rectangle
            %   width_x
            %       Desc: scalar double, width of rectangle
            %   height_y
            %       Desc: scalar double, height of rectangle
            %   index
            %       Desc: index of refraction
            
            x = obj.x_coords;
            y = obj.y_coords;
            
            % fill in the rect
            obj.N( y > (min_y - obj.dy/2) & y <= (min_y + height_y - obj.dy/2), ...
                    x > (min_x - obj.dx/2) & x <= (min_x + width_x - obj.dx/2) ) = index; 
            
        end     % end function addRect()
        
        
        function obj = twoLevelBuilder( obj, wgs_min_y, wgs_thick, wgs_indx, ...
                                             wgs_duty_cycles, wgs_offsets )
            % Draws the two level waveguide grating
            % 
            % Inputs:
            %   wg_min_y
            %       type: double, array
            %       desc: 1x2 array, bottom y pos of each level
            %   wg_thick
            %       type: double, array
            %       desc: 1x2 array, thickness of each level
            %   wg_indx
            %       type: double, array
            %       desc: 1x2 array, index of each level
            %   wgs_duty_cycles
            %       type: double, array
            %       desc: 1x2 array, duty cycle of each level, defined as
            %             waveguide length/unit cell period
            %   wgs_offsets
            %       type: double, array
            %       desc: 1x2 array, offset of each grating tooth from the
            %             left edge of the unit cell
            %
            % Sets these properties:
            %   obj.wg_min_y
            %   obj.wg_max_y
            
            
            % grab period
            a = obj.domain_size(2);
            
            % add first grating tooth
            tooth1_length   = a*wgs_duty_cycles(1);
            obj             = obj.addRect( wgs_offsets(1), wgs_min_y(1), tooth1_length, wgs_thick(1), wgs_indx(1) );
            if ( wgs_offsets(1) + tooth1_length > a )
                % If the tooth size extends past the domain, wrap it back
                % around from the start
                
                tooth1_length_remainder = wgs_offsets(1) + tooth1_length - a;
                obj                     = obj.addRect( 0, wgs_min_y(1), tooth1_length_remainder, wgs_thick(1), wgs_indx(1) );
                
            end
            
            % add 2nd grating tooth
            tooth2_length   = a*wgs_duty_cycles(2);
            obj             = obj.addRect( wgs_offsets(2), wgs_min_y(2), tooth2_length, wgs_thick(2), wgs_indx(2) );
            if ( wgs_offsets(2) + tooth2_length > a )
                % If the tooth size extends past the domain, wrap it back
                % around from the start
                
                tooth2_length_remainder = wgs_offsets(2) + tooth2_length - a;
                obj                     = obj.addRect( 0, wgs_min_y(2), tooth2_length_remainder, wgs_thick(2), wgs_indx(2) );
                
            end
            
            % set top and bottom bounds for guided region
            obj.wg_min_y                = min( wgs_min_y );     % bottom position of wg
            [ obj.wg_max_y, indx_max ]  = max( wgs_min_y );     % top position of wg
            obj.wg_max_y                = obj.wg_max_y + wgs_thick( indx_max );
            
        end
        
        
        function obj = runSimulation( obj, num_modes, BC, pml_options, guessk, OPTS )
            % Runs new mode solver
            %
            % Description:
            %   Runs complex-k mode solver. Stores the mode with the most
            %   guided power. Calculates up/down power, directivity,
            %   angle of maximum radiation, and scattering strength.
            %
            % Inputs:
            %   num_modes
            %       type: integer
            %       desc: # of modes (max) to simulate
            %   BC
            %       type: integer
            %       desc: 0 for PEC, 1 for PMC
            %   pml_options
            %       type: array, double
            %       desc: 1x4 Array with the following elements:
            %               PML_options(1): PML in y direction (yes=1 or no=0)
            %               PML_options(2): length of PML layer in nm
            %               PML_options(3): strength of PML in the complex plane
            %               PML_options(4): PML polynomial order (1, 2, 3...)
            %   guessk
            %       type: scalar, double (can be complex)
            %       desc: guess k value. Works best when closest to desired
            %             mode. In units rad/'units'
            %   OPTS
            %       type: struct
            %       desc: optional options with the following fields
            %           'mode_to_overlap'
            %               type: matrix, double
            %               desc: mode to overlap
            %
            % Sets these properties:
            %   obj.k
            %       units rad/'units'
            %   obj.Phi
            %       dimensions y vs. x (x is dir. of propagation)
            %   obj.E_z
            %       field repeated, i can't remember why or if i use this
            %   obj.directivity
            %       up/down power ratio
            %   calls obj.calc_output_angle()
            %       this is a function, which probably figures out the
            %       output angle
            %   calls obj.calc_scattering_strength()
            %

            % default OPTS
            if nargin < 6
                OPTS = struct();
            end
            
            % spatial variables, in units nm
            nm          = 1e9;
            a           = obj.domain_size(2) * obj.units.scale * nm;
            lambda_nm   = obj.lambda * obj.units.scale * nm;
            dx_nm       = obj.dx * obj.units.scale * nm;
            guessk_nm   = guessk / ( obj.units.scale * nm );                % units rad/nm

            % store options
            obj.sim_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options, 'OPTS', OPTS );

            % set guessk if not entered
            if nargin < 5
                guessk = pi/(2*a);
            end
            
            % run solver
            k0_nm           = 2*pi/lambda_nm;
            [Phi_all, k_nm] = complexk_mode_solver_2D_PML( obj.N, ...
                                                       dx_nm, ...
                                                       k0_nm, ...
                                                       num_modes, ...
                                                       guessk_nm, ...
                                                       BC, ...
                                                       pml_options );
                                                   
            % re-scale k to units 'units'
            k_vs_mode = k_nm * nm * obj.units.scale;
            
            % save k and phi vs mode #
            obj.k_vs_mode   = k_vs_mode;
            obj.Phi_vs_mode = Phi_all;
            
%             % DEBUG store temporary copies of k and phi_all b4 removing and
%             % sorting
%             obj.debug.k_all_before_fix_backwards     = k;
%             obj.debug.phi_all_before_fix_backwards   = Phi_all;

            
%             % sort on guided power
%             guided_power  = zeros( size(k) );
%             total_power   = zeros( size(k) );
% 
%             % check to see if waveguide boundaries have been set yet
%             if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
%                 error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
%                         ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
%             end
%             
%             y_bot   = obj.wg_min_y;
%             y_top   = obj.wg_max_y;
%             y       = obj.y_coords;
%             
%             
%             % Calculated guided power and fix backwards modes if option is
%             % selected
%             for ii = 1:length(k)
                % for each mode
% 
%                 % Fix backwards propagating modes
%                 if isfield( OPTS, 'fix_neg_k' ) && OPTS.fix_neg_k == true
% 
%                     % check if k is backwards propagating
%                     % if it is, then the mode must be resimulated with -k
%                     if imag( k(ii) ) < 0
% 
%                         fprintf([ '\nMode found is backwards propagating (positive real k, negative imag k).\n', ...
%                                   'Re-running solver with flipped k\n' ]);
%                         fprintf('Current k = %f + i * %f\n', real( k(ii) ), imag( k(ii) ) );
% 
%                         % re-run solver
%                         [Phi, k_nm] = complexk_mode_solver_2D_PML( obj.N, ...
%                                                                    dx, ...
%                                                                    k0, ...
%                                                                    1, ...
%                                                                    -( k(ii) * obj.units.scale * nm ), ...
%                                                                    BC, ...
%                                                                    pml_options );
% 
%                         % re-scale k
%                         k(ii) = k_nm * nm * obj.units.scale;
% 
%                         % now check, if imag(k) is still < 0, kill the mode
%                         % from consideration, meaning set the field = 0
%                         if imag( k(ii) ) < 0
%                             % set field = 0
%                             fprintf('Mode %i has imag(k) < 0, killing this mode\n', ii );
%                             Phi_all( :, :, ii ) = zeros(size(Phi));
%                         else
%                             % overwrite Phi
%                             fprintf('New k = %f + i * %f\n', real( k(ii) ), imag( k(ii) ) );
%                             Phi_all( :, :, ii ) = Phi;
%                         end
% 
%                     end     % end if imag( k(ii) ) < 0
% 
%                 end     % end if isfield
%                     
%                 % grab guided portion of field
%                 phi_guided  = Phi_all( y >= y_bot & y <= y_top, :, ii );
%                 cur_phi     = Phi_all( :, :, ii );
% 
%                 % sum area of field
%                 guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
%                 total_power(ii)     = sum( abs( cur_phi(:).^2 ) );
% 
%             end     % end for
  
            
            
            
%             % DEBUG storing the guided power
%             obj.debug.guided_power          = guided_power;
%             obj.debug.guided_power_ratio    = guided_power./total_power;
%             
%             % keep mode with LEAST unguided power OR MOST guided power
%             [~, indx_k] = max( abs(guided_power./total_power) );        % most guided
%             k           = k(indx_k);
%             Phi         = Phi_all(:,:,indx_k);
%             
% 
%             % save wavevectors and field
%             obj.k               = k;
%             obj.Phi             = Phi;
%             obj.chosen_mode_num = indx_k;
            
            % pick which mode to keep
            if isfield( OPTS, 'mode_to_overlap' )
                % pick mode with best overlap
                obj = obj.choose_mode( OPTS.mode_to_overlap );
            else
                % pick mode guided mode
                obj = obj.choose_mode(); 
            end
            
            % stitch together full e field, with the request number of
            % periods
            [obj, E_z]  = obj.stitch_E_field( obj.Phi, obj.k, obj.numcells );
            obj.E_z     = E_z;
            
            % for mode overlapping, stitch together single unit cell of E
            % field, using only real(k)
            [obj, E_z_for_overlap]  = obj.stitch_E_field( obj.Phi, real(obj.k), 1 );
            obj.E_z_for_overlap     = E_z_for_overlap;
            
            % pick slices of field to compute directivity, angle, etc.
            h_pml_d = round( pml_options(2)/obj.dy );                       % size of pml in discretizations
            y_up    = size( obj.E_z, 1 ) - h_pml_d - 1;
            y_down  = h_pml_d+2;
            
            % calculate output angle
            obj = obj.calc_output_angle( y_up, y_down );
            
            % calculate up/down directivity
            obj             = obj.calc_radiated_power();
            obj.directivity = obj.P_rad_up/obj.P_rad_down;
                    
            % calculate power scattering strength
            obj = obj.calc_scattering_strength();
            
            
        end     % end function runSimulation()
        
        
        function obj = choose_mode( obj, mode_to_overlap )
            % Function that chooses which mode becomes the accepted mode
            %
            % Inputs:
            %   mode_to_overlap
            %       type: matrix, double
            %       desc: OPTIONAL. If passed in as argument, this function
            %             will choose the mode with the closest overlap.
            %             Otherwise, the function will choose the mode with
            %             the mode guided power
            %
            % Sets these properties:
            %   obj.k
            %       prop. k of chosen mode
            %   obj.Phi
            %       field envelope of chosen mode
            %   obj.chosen_mode_num
            %       index of chosen mode (chosen mode k = obj.k_vs_mode(
            %       obj.chosen_mode_num))
            
            if nargin < 2
                % sort on guided power
                
                guided_power  = zeros( size(obj.k_vs_mode) );
                total_power   = zeros( size(obj.k_vs_mode) );

                % check to see if waveguide boundaries have been set yet
                if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
                    error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
                            ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
                end

                y_bot   = obj.wg_min_y;
                y_top   = obj.wg_max_y;
                y       = obj.y_coords;

                % Calculated guided power and fix backwards modes if option is
                % selected
                for ii = 1:length(obj.k_vs_mode)
                    % for each mode

                    % grab guided portion of field
                    phi_guided  = obj.Phi_vs_mode( y >= y_bot & y <= y_top, :, ii );
                    cur_phi     = obj.Phi_vs_mode( :, :, ii );

                    % sum area of field
                    guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
                    total_power(ii)     = sum( abs( cur_phi(:).^2 ) );

                end     % end for

                % DEBUG storing the guided power
                obj.debug.guided_power          = guided_power;
                obj.debug.guided_power_ratio    = guided_power./total_power;

                % keep mode with LEAST unguided power OR MOST guided power
                [~, indx_k]         = max( abs(guided_power./total_power) );        % most guided
                obj.k               = obj.k_vs_mode(indx_k);
                obj.Phi             = obj.Phi_vs_mode(:,:,indx_k);
                obj.chosen_mode_num = indx_k;
                
            else
                % sort on mode overlap
                
                % run overlaps
                [obj, max_overlaps] = obj.calc_mode_overlaps( mode_to_overlap );
                
                % keep mode with highest overlap
                [~, indx_k]         = max( max_overlaps );        
                obj.k               = obj.k_vs_mode(indx_k);
                obj.Phi             = obj.Phi_vs_mode(:,:,indx_k);
                obj.chosen_mode_num = indx_k;
                
            end
            
        end     % end function choose_mode()
        
        
        function [obj, max_overlaps] = calc_mode_overlaps( obj, mode_to_overlap )
            % Takes in a mode and overlaps it with the modes that were
            % solved in this object
            % Overlap is performed using a cross correlation
            %
            % Inputs:
            %   mode_to_overlap
            %       type: matrix, double
            %       desc: field to overlap with. Remove the exponential
            %             term from it
            %
            % Outputs:
            %   max_overlaps
            %       type: vector, double
            %       desc: Max mode overlap vs. mode #     
            
            % number of modes
            n_modes = length( obj.k_vs_mode );
            
            % size of mode to overlap
            [ ny_mode_to_overlap, nx_mode_to_overlap ] = size( mode_to_overlap );

            % save max overlaps
            max_overlaps = zeros( n_modes, 1 );
            
            for mode_num = 1:n_modes
                % for each of this object's saved modes
                
                % grab size of this object's mode
                [ ny_this_mode, nx_this_mode ]  = size( obj.Phi_vs_mode(:,:,mode_num) );
            
                % grab this mode, including phase
                this_mode   = obj.Phi_vs_mode(:,:,mode_num) .* repmat( exp( 1i * obj.x_coords * real(obj.k_vs_mode( mode_num ) )), ny_this_mode, 1 );
            
                % zero pad
%                 total_ny            = ny_this_mode + ny_mode_to_overlap;
                mode_to_overlap_pad = padarray( mode_to_overlap, [ floor(ny_this_mode/2), floor(nx_this_mode/2) ], 0, 'pre' );
                mode_to_overlap_pad = padarray( mode_to_overlap_pad, [ ceil(ny_this_mode/2), ceil(nx_this_mode/2) ], 0, 'post' );
                this_mode_pad       = padarray( this_mode, [ floor(ny_mode_to_overlap/2), floor(nx_mode_to_overlap/2) ], 0, 'pre' );
                this_mode_pad       = padarray( this_mode_pad, [ ceil(ny_mode_to_overlap/2), ceil(nx_mode_to_overlap/2) ], 0, 'post' );
            
                % x-correlate
                mode_xcorr = ifftshift( ifft2( conj( fft2( mode_to_overlap_pad ) ) .* fft2( this_mode_pad ) ) );
                
                % save max overlap
                max_overlaps( mode_num ) = max( abs( mode_xcorr(:) ) );

            end

        end     % end function calc_mode_overlaps()
        
        

        
        
        function obj = calc_radiated_power(obj)
            % Calculates the radiated power in both the upwards and
            % downwards directions
            % Only run AFTER mode solver has been run, since this function
            % depends on E_z and k
            %
            % Inputs:
            %   y_up    - location of top slice of E field to take
            %   y_down  - location of bot slice of E field to take
            %       y_up and y_down are deprecated.
            %
            % Uses these object properties:
            %   E_z
            %   k
            %
            % Sets these object properties:
            %   P_rad_down      - saves power radiated down
            %   P_rad_up        - saves power radiated up
            
%             % define constants
%             mu0     = 4*pi*1e-7;                % units of H/m
%             mu0     = mu0 * obj.units.scale;    % units of H/(units)
%             c       = 3e8;                      % units of m/s
%             c       = c/obj.units.scale;        % units of (units)/s
%             omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
%             
%            
%             % Stich field and phase together
%             phase_onecell   = repmat( exp( 1i * obj.k * obj.x_coords ), size( obj.Phi, 1 ), 1 );
%             Ez_onecell      = obj.Phi .* phase_onecell;
%            
%             
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
%             
%             % save to obj
%             obj.Sx              = Sx;
%             obj.P_per_x_slice 	= P_per_x_slice;
%             
% 
%             % Calculate power radiated up
%             
%             
%             % Calculate and plot power propagating upwards vs. y
%             
%             % calculate H_x
%             % dimensions H_x vs. y vs x
%             H_x  = 1/(1i*omega0*mu0) .* ( Ez_onecell( 3:end,:) - Ez_onecell( 1:end-2,:) )/(2*obj.dy);
% 
%             % power flowing across y
%             Sy                  = real( conj( H_x ) .* Ez_onecell( 2:end-1,: ) );       % Sy = real( Ez Hx* )
%             P_per_y_slice       = sum( Sy, 2 )*obj.dx;                                  % using Sy compmonent
            
            % calculate power distributions for the chosen mode
            [obj, Sx, Sy, P_per_x_slice, P_per_y_slice] = obj.calc_power_distribution_onecell();
            
            
            [~, indx_max_up]    = max( P_per_y_slice );
            [~, indx_max_down]  = min( P_per_y_slice );
            Sy_up               = Sy( indx_max_up, : );                                 % dimensions Sy vs. x
            Sy_down             = Sy( indx_max_down, : );
            
            % save to obj
            obj.Sx              = Sx;
            obj.P_per_x_slice 	= P_per_x_slice;
            obj.Sy              = Sy;
            obj.P_per_y_slice   = P_per_y_slice;
            
        
%             % new version, take the flux including both Sx and Sy
%             % components
%             P_rad_up    = sum( sqrt( abs( Sy_up(2:end-1) ).^2 + abs( Sx(y_up,:) ).^2 ) ) * obj.dx;
%             P_rad_down  = sum( sqrt( abs( Sy_down(2:end-1) ).^2 + abs( Sx(y_down,:) ).^2 ) ) * obj.dx;

            % new new version, take only Sy but multiply by ./cos(theta)
            P_rad_up    = sum( abs( Sy_up(:) ) ) * obj.dx ./ cos( obj.max_angle_up * pi/180 );
            P_rad_down  = sum( abs( Sy_down(:) ) ) * obj.dx ./ cos( obj.max_angle_down * pi/180 );
            
            % save to object
            obj.P_rad_down    = P_rad_down;
            obj.P_rad_up      = P_rad_up;
            
            
            % DEBUG calculate angle from Sx and Sy
% %             obj.debug.tan_sx_sy_up = (180/pi) * atan( Sx(y_up,:)./ Sy_up(2:end-1) );
        
            
        end     % end function calc_radiated_power()
        
        
        function obj = calc_output_angle(obj, y_up, y_down)
            % Calculates output angle of radiation
            % To be run only AFTER mode solver has been run
            %
            % Inputs:
            %   y_up    - location of top slice of E field to take
            %   y_down  - location of bot slice of E field to take
            
            % load stuff
%             E_z = obj.E_z;
            f0  = 1/obj.lambda;     % spatial freq
            
            % generate field, without the exponential growth or decay
            numcells        = 2000;                                                  % i'm purposefully setting my own number of cells here
            nx              = round( numcells*obj.domain_size(2)/obj.dx );
            x_coords_all    = 0 : obj.dx : ( (nx-1) * obj.dx );
            phase_all       = repmat( exp( 1i*real(obj.k)*x_coords_all ), size(obj.Phi,1), 1 );
            E_z             = repmat( obj.Phi, 1, numcells ).*phase_all;
            
            % grab slices
            E_z_up          = E_z( y_up, : );
            E_z_down        = E_z( y_down, : );
            
            % do some zero padding and repeating
            % currently zero padding is turned off.
%             pad_factor  = 0;
%             E_z_up      = [ zeros( 1, pad_factor*length(E_z_up) ), E_z_up, zeros( 1, pad_factor*length(E_z_up) ) ];
%             E_z_down    = [ zeros( 1, pad_factor*length(E_z_down) ), E_z_down, zeros( 1, pad_factor*length(E_z_down) ) ];
            
            % index of refraction of top and bottom cladding
            n_top = obj.N( y_up, 1 );
            n_bot = obj.N( y_down, 1 );
                 
            % create space and frequency vectors
            Nx      = length( E_z_up );
            fx_vec  = (-Nx/2:(Nx/2-1) ).*(1/(obj.dx*Nx));      % spatial freq vector

            % vectors that relate fz to angle, in deg
            angle_vec_bot               = (180/pi) * asin( fx_vec./(n_bot*f0) );    % draw a circ
            real_a_indices_bot          = abs(imag(angle_vec_bot)) < 1e-6 ;         % use this to index the range of allowed free space modes
            angle_vec_top               = (180/pi) * asin( fx_vec./(n_top*f0) );    % draw a circ
            real_a_indices_top          = abs(imag(angle_vec_top)) < 1e-6 ;         % use this to index the range of allowed free space modes

            % take FFT
            E_z_fx_up       = fftshift( fft( ifftshift( E_z_up ) ) );
            E_z_fx_down     = fftshift( fft( ifftshift( E_z_down ) ) );
            
%             % DEBUG plot E_fz
%             % first normalize E fields
%             E_z_fx_down_norm    = E_z_fx_down./max(abs(E_z_fx_down(:)));
%             E_z_fx_up_norm      = E_z_fx_up./max(abs(E_z_fx_up(:)));
%             
%             % DEBUG plot E_fz down
%             figure;
%             plot( angle_vec_bot(real_a_indices_bot), abs(E_z_fx_down_norm(real_a_indices_bot)) );
%             xlabel('Rad. angle (deg)'); ylabel('E_z(f_x)');
%             title('DEBUG abs(E(f_x)) down, within radiation window');
%             makeFigureNice();
%             
%             % DEBUG plot E_fz down
%             fy_down     = real(sqrt( (f0*n_bot)^2 - fx_vec.^2 ));     % fy component
%             fy_down     = fy_down./max(abs(fy_down));
%             figure;
%             plot( fx_vec, abs(E_z_fx_down_norm) ); hold on;
%             plot( fx_vec, fy_down, '--' );
%             xlabel('f_x (1/\lambda)'); ylabel('E_z(f_x)');
%             legend('E_z', 'rad. circle');
%             title('DEBUG abs(E(f_x)) down, ALL');
%             makeFigureNice();
%             
%             % DEBUG plot E_fz up
%             figure;
%             plot( angle_vec_top(real_a_indices_top), abs(E_z_fx_up_norm(real_a_indices_top)) );
%             xlabel('Rad. angle (deg)'); ylabel('E_z(f_x)');
%             title('DEBUG abs(E(f_x)) up, within radiation window');
%             makeFigureNice();
%             
%             % DEBUG plot E_fz up
%             fy_up     = real(sqrt( (f0*n_top)^2 - fx_vec.^2 ));     % fy component
%             fy_up     = fy_up./max(abs(fy_up));
%             figure;
%             plot( fx_vec, abs(E_z_fx_up_norm) ); hold on;
%             plot( fx_vec, fy_up, '--' );
%             xlabel('f_x (1/\lambda)'); ylabel('E_z(f_x)');
%             legend('E_z', 'rad. circle');
%             title('DEBUG abs(E(f_x)) up, ALL');
%             makeFigureNice();
            
            % calc and save max output angle
            % going up
            E_z_fx_up           = E_z_fx_up(real_a_indices_top);
            angle_vec_top       = angle_vec_top(real_a_indices_top);
            [~, indx_max_up]    = max( abs( E_z_fx_up ) );
            obj.max_angle_up    = angle_vec_top(indx_max_up);
            % going down
            E_z_fx_down         = E_z_fx_down(real_a_indices_bot);
            angle_vec_bot       = angle_vec_bot(real_a_indices_bot);
            [~, indx_max_down]  = max( abs( E_z_fx_down ) );
            obj.max_angle_down  = angle_vec_bot(indx_max_down);      
            
        end     % end function calc_output_angle()
        
        
        function obj = calc_scattering_strength(obj)
            % GOING TO REINSTATE THIS FUNCTION
            % Calculates scattering strength in terms of decay rate alpha
            % To be run only AFTER mode solver AND calc_radiated_power have
            % been run
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            mu0     = mu0 * obj.units.scale;    % units of H/(units)
            c       = 3e8;                      % units of m/s
            c       = c/obj.units.scale;        % units of (units)/s
            omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
            
            % grab field
            E_z = obj.E_z;
            
            % H x in
            % gonna calculate input power at grid point 2
            H_x_in1 = (-1i/(omega0*mu0)) * ( E_z( 3:2:end, 2 ) - E_z( 1:2:end-2, 2 ) )/(2*obj.dy);   % dy, term staggered 1
            H_x_in2 = (-1i/(omega0*mu0)) * ( E_z( 4:2:end, 2 ) - E_z( 2:2:end-2, 2 ) )/(2*obj.dy);   % dy, term staggered 2
            % interleave two arrays
            H_x_in                                  = zeros( size( [H_x_in1, H_x_in2] ) );
            H_x_in( ( 1:length(H_x_in1) )*2 - 1)    = H_x_in1;
            H_x_in( ( 1:length(H_x_in2) )*2 )       = H_x_in2;
            
            % H y in (on the same grid as Hx in)
%             H_y_in = (1i/(omega0*mu0)) * ( E_z( 2:end-1, 3 ) - E_z( 2:end-1, 1 ) )/(2*obj.dx);   % dy, term staggered 1
            
            % NEW: H y in , on SECOND step, using entire E_z slice
            H_y_in = (1i/(omega0*mu0)) * ( E_z( :, 3 ) - E_z( :, 1 ) )/(2*obj.dx);   % dy, term staggered 1
            
            % power in (at left edge)
            Sx_in   = real( -1 * conj(H_y_in(:)) .* E_z( :, 2 ) );              % Sx = real( -Ez Hy* )
            P_in    = sum( Sx_in(:) )*obj.dy;

            
            % DEBUG save Sx_in
            obj.debug.Sx_in = Sx_in;

            % save input power
            obj.P_in = P_in;
            
            % calculate alpha efficiency (in 1/units)
            period          = obj.domain_size(2);
            obj.alpha_up    = ( -1/( 2*period ) ) * log( 1 - obj.P_rad_up/P_in );
            obj.alpha_down  = ( -1/( 2*period ) ) * log( 1 - obj.P_rad_down/P_in );


            % DEBUG calculate efficiency from one cell
            
                     
%             % DEBUG let's calculate, numerically, the loss coefficient from
%             % input to output and compare with imag(k)
%             % H x out, at grid point end-1
%             H_x_out1 = (-1i/(omega0*mu0)) * ( E_z( 3:2:end, end-1 ) - E_z( 1:2:end-2, end-1 ) )/(2*obj.dy);   % dy, term staggered 1
%             H_x_out2 = (-1i/(omega0*mu0)) * ( E_z( 4:2:end, end-1 ) - E_z( 2:2:end-2, end-1 ) )/(2*obj.dy);   % dy, term staggered 2
%             % interleave two arrays
%             H_x_out                                 = zeros( size( [H_x_out1, H_x_out2] ) );
%             H_x_out( ( 1:length(H_x_out1) )*2 - 1)  = H_x_out1;
%             H_x_out( ( 1:length(H_x_out2) )*2 )     = H_x_out2;
%             
%             % H y out (on the same grid as Hx out)
%             H_y_out = (1i/(omega0*mu0)) * ( E_z( 2:end-1, end ) - E_z( 2:end-1, end-2 ) )/(2*obj.dx);   % dy, term staggered 1
%             
%             % power in (at left edge)
%             Sy_out = real( conj(H_x_out(:)) .* E_z( 2:end-1, end-1 ) );         % Sy = real( Ez Hx* )
%             Sx_out = real( -1 * conj(H_y_out(:)) .* E_z( 2:end-1, end-1 ) );    % Sx = real( -Ez Hy* )
%             P_out  = sum( sqrt( Sy_out.^2 + Sx_out.^2 ) )*obj.dy;                 % using both Sx and Sy compmonents
%             
%             % numerical calculation of loss coeff
%             alpha_tot = ( -1/( 2*obj.numcells*period ) ) * log( P_out/P_in );
            
        end     % end function calc_scattering_strength()
        
        
        function obj = calc_power_distribution(obj)
            % Mostly for debugging
            % plots x y dependence of power
            %
            % also saving poynting vector components Sx and Sy to
            % obj.debug.Sx and obj.debug.Sy
            % Sx and Sy have same dimensions as E field - y vs. x
            % except Sx is only defined from x_coords(2:end-1)
            % and Sy is only defined from y_coords(2:end-1)
            % and saving power per slice as:
            % obj.debug.P_per_y_slice, obj.debug.P_per_x_slice
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            mu0     = mu0 * obj.units.scale;    % units of H/(units)
            c       = 3e8;                      % units of m/s
            c       = c/obj.units.scale;        % units of (units)/s
            omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
            
            % grab field
            E_z             = obj.E_z;
            x_coords_all    = 0 : obj.dx : obj.numcells*obj.domain_size(2)-obj.dx;
            obj.debug.x_coords_all = x_coords_all;
            
            % Calc and plot propagating power vs. x
            
            % H y in , on second to end-1 steps, using entire E_z
            % dimensions are y (transverse) vs. x
            H_y_in = (1i/(omega0*mu0)) * ( E_z( :, 3:end ) - E_z( :, 1:end-2 ) )/(2*obj.dx);   % dy, term staggered 1
            
            % poynting vector in, dimensions are y vs x
            Sx = real( -1 * conj( H_y_in ) .* E_z( :, 2:end-1 ) );              % Sx = real( -Ez Hy* )

            % power per x slice
            P_per_x_slice = sum( Sx, 1 )*obj.dy;
            
            % lets plot power traveling in x per slice of x
            figure;
            plot( x_coords_all(2:end-1), P_per_x_slice );
            xlabel('x'); ylabel('total propagating power in x');
            title('Power propagating along x per slice');
            makeFigureNice();
           
            
            % Calculate and plot power propagating upwards vs. y
            
            % calculate H_x
            % dimensions H_x vs. y vs x
            H_x  = 1/(1i*omega0*mu0) .* ( E_z( 3:end,:) - E_z( 1:end-2,:) )/(2*obj.dy);

            % power flowing across y
            Sy              = real( conj( H_x ) .* E_z( 2:end-1,: ) );      % Sy = real( Ez Hx* )
            P_per_y_slice   = sum( Sy, 2 )*obj.dx;                          % using Sy compmonent
            
            % plot power traveling in y per slice of y
            figure;
            plot( obj.y_coords(2:end-1), P_per_y_slice );
            xlabel('y'); ylabel('total propagating power in y');
            title('Power propagating along y per slice');
            makeFigureNice();
            
            % DEBUG Save Sx and Sy, power per slice
            obj.debug.Sx            = Sx;
            obj.debug.Sy            = Sy;
            obj.debug.P_per_y_slice = P_per_y_slice;
            obj.debug.P_per_x_slice = P_per_x_slice;
            
        end     % end function calc_power_distribution()
        
        
        function [obj, Sx, Sy, P_per_x_slice, P_per_y_slice] = ...
                        calc_power_distribution_onecell(obj, mode_num)
            % Calculates Sx and Sy, power distribution for a single unit
            % cell
            %
            % Pre-requisite for running calc_radiated_power
            %
            % inputs:
            %   mode_num
            %       type: scalar, int
            %       desc: OPTIONAL, mode number to calc Sx and Sy for
            %
            % outputs:
            %   Sx
            %       type: matrix, double
            %       desc: x component (in dir. of waveguide prop) of poynting vector
            %             dimensions y vs. x(2:end-1)
            %   Sy
            %       type: matrix, double
            %       desc: y component (transverse dir, pos = up) of poynting vector
            %             dimensions y(2:end-1) vs. x
            %   P_per_x_slice
            %       type: vector, double
            %       desc: integral of Sx aka power propagating through each
            %             x slice of the domain
            %             dimensions power vs. x(2:end-1)
            %   P_per_y_slice
            %       type: vector, double
            %       desc: integral of Sy aka power propagating through each
            %             y slice of the domain
            %             dimensions power vs. y(2:end-1)
            
            if nargin < 2
                choose_mode = false;
            else
                choose_mode = true;
            end
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            mu0     = mu0 * obj.units.scale;    % units of H/(units)
            c       = 3e8;                      % units of m/s
            c       = c/obj.units.scale;        % units of (units)/s
            omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
            
            % Choose field to use
            if ~choose_mode
                % use the chosen "most guided" mode

                % Stich field and phase together
                phase_onecell   = repmat( exp( 1i * obj.k * obj.x_coords ), size( obj.Phi, 1 ), 1 );
                Ez_onecell      = obj.Phi .* phase_onecell;
            else
                % choose which mode to use
                
                k_chosen    = obj.k_vs_mode( mode_num );
                phi_chosen  = obj.Phi_vs_mode( mode_num );
                % Stich field and phase together
                phase_onecell   = repmat( exp( 1i * k_chosen * obj.x_coords ), size( phi_chosen, 1 ), 1 );
                Ez_onecell      = phi_chosen .* phase_onecell;
                
            end
       
            % Calc propagating power vs. x
            
            % H y in , on second to end-1 steps, using entire E_z
            % dimensions are y (transverse) vs. x
            H_y_in = (1i/(omega0*mu0)) * ( Ez_onecell( :, 3:end ) - Ez_onecell( :, 1:end-2 ) )/(2*obj.dx);   % dx, term staggered 1
            
            % poynting vector in, dimensions are y vs x(2:end-1)
            Sx = real( -1 * conj( H_y_in ) .* Ez_onecell( :, 2:end-1 ) );              % Sx = real( -Ez Hy* )

            % power per x slice
            P_per_x_slice = sum( Sx, 1 )*obj.dy;
            

            % Calculate power radiated up
            
            % calculate H_x
            % dimensions H_x vs. y vs x
            H_x  = 1/(1i*omega0*mu0) .* ( Ez_onecell( 3:end,:) - Ez_onecell( 1:end-2,:) )/(2*obj.dy);

            % power flowing across y
            Sy                  = real( conj( H_x ) .* Ez_onecell( 2:end-1,: ) );       % Sy = real( Ez Hx* )
            P_per_y_slice       = sum( Sy, 2 )*obj.dx;                                  % using Sy compmonent

            
        end     % end function calc_power_distribution_onecell()
        
        
        function [obj, E_z] = stitch_E_field( obj, Phi, k, num_cells )
            % Stitches the E field together from the phase and envelope
            %
            % Inputs:
            %   Phi
            %       type: matrix, double
            %       desc: field envelope
            %   k
            %       type: scalar, double
            %       desc: propagation constant
            %   num_cells
            %       type: scalar, int
            %       desc: number of cells to repeat
            
            % stitch together e field, including the phase
            nx              = round( num_cells*obj.domain_size(2)/obj.dx );
            x_coords_all    = 0 : obj.dx : ( (nx-1) * obj.dx );
            phase_all       = repmat( exp( 1i*k*x_coords_all ), size(Phi,1), 1 );
            E_z             = repmat( Phi, 1, num_cells ).*phase_all;
            
        end     % end function stitch_E_field()
        
        
        function obj = shift_index_circ( obj, shift_length_x )
            % Circularly shifts the index by the shift length
            % Shift length can be positive or negative
            % Ideally the shift length is an integer multiple of obj.dx
            %
            % Inputs:
            %   shift_length_x
            %       type: double, scalar
            %       desc: length to shift by, can be either positive or
            %             negative value
            
            % calc number of units to shift by
            nx      = round( shift_length_x/obj.dx );
            obj.N   = circshift( obj.N, nx, 2 );
            
        end     % end function shift_index_circ()
        
        
        function plotIndex(obj)
            % Plots the index distribution
            
            figure;
            imagesc( obj.x_coords, obj.y_coords, obj.N );
            colorbar;
            set( gca, 'YDir', 'normal' );
            xlabel(['x (', obj.units.name, ')']);
            ylabel(['y (', obj.units.name, ')']);
            title('Plot of index distribution');
            
        end     % end function plotIndex()
        
        
        function plotEz(obj)
            % plots electric field, this could probably be part of a gui
            
            numcells = obj.numcells;
            E_z      = obj.E_z;
            
            % repmat the index
            N = repmat( obj.N, 1, numcells );
            
            % define x coordinates
            x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
            
            
            % plot the field, real
            figure;
            imagesc( x_coords_all, obj.y_coords, real(E_z) );
            clim( [ -max( abs( real(E_z(:)) ) ), max( abs( real(E_z(:)) ) ) ] );                % set color limits
            xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
            set( gca, 'YDir', 'normal' );
            title( sprintf('Real E_z, %i periods', numcells) );
            
            % plot the field, abs
            figure;
            imagesc( x_coords_all, obj.y_coords, abs(E_z) );
            xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
            set( gca, 'YDir', 'normal' );
            title( sprintf('Amplitude E_z, %i periods', numcells) );
            

            
        end     % end function plotEz()
        
        
        function plotEz_w_edges(obj, OPTS)
            % plots electric field, this could probably be part of a gui
            % with grating geometry edges overlaid
            %
            % USE THE GUI INSTEAD
            % obj.plot_E_field_gui()
            %
            % inputs:
            %   OPTS
            %       type: struct
            %       desc: struct with these fields:
            %           plots
            %               type: string
            %               desc: either 'real', 'amp', or 'both'
            %           mode_num
            %               type: scalar, int
            %               desc: Which mode to plot
            %                     Mainly for debugging
            %                     Optional, if not set, then the plotted
            %                     mode is the chosen "most guided" mode
            
            % default to both
            if nargin < 2
                OPTS = struct( 'plots', 'both' );
            end
            if ~isfield( OPTS, 'plots' )
                OPTS.plots = 'both';
            end
            
            numcells = obj.numcells;
            E_z      = obj.E_z;
            
            % repmat the index
            N = repmat( obj.N, 1, numcells );
            
            % define x coordinates
            x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
            
            
            if ~isfield( OPTS, 'mode_num' )
                % Plot the chosen "most guided" mode
                field_to_plot = E_z;
            else
                % plot the selected mode #
                
                % stitch together e field, including the phase
                phi_chosen_mode = obj.Phi_vs_mode( :, :, OPTS.mode_num );
                k_chosen_mode   = obj.k_vs_mode( OPTS.mode_num );
                phase_all       = repmat( exp( 1i*k_chosen_mode*x_coords_all ), size(phi_chosen_mode,1), 1 );
                field_to_plot   = repmat( obj.Phi_vs_mode(:,:,OPTS.mode_num), 1, obj.numcells ).*phase_all;
                
            end
            
            % plot the field, real
            if strcmp( OPTS.plots, 'real' ) || strcmp( OPTS.plots, 'both' ) 
                figure;
                imagesc( x_coords_all, obj.y_coords, real(field_to_plot) ); hold on;
                xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
                set( gca, 'YDir', 'normal' );
                title( sprintf('Real E_z, %i periods', numcells) );
                % overplot the index with transparency
                % first things first, edge detection:
                filt    = 'roberts';            % filter
                N_edges = edge( N, filt );
                % overplot the detected edges
                h       = imagesc( x_coords_all, obj.y_coords, repmat( ~N_edges, 1, 1, 3 ) );
                set( h, 'AlphaData', N_edges );
                set( gca, 'YDir', 'normal' );
            end

            % plot the field, abs
            if strcmp( OPTS.plots, 'amp' ) || strcmp( OPTS.plots, 'both' ) 
                figure;
                imagesc( x_coords_all, obj.y_coords, abs(field_to_plot) ); hold on;
                xlabel('x'); ylabel('y'); colormap('redbluehilight'); colorbar;
                set( gca, 'YDir', 'normal' );
                title( sprintf('Amplitude E_z, %i periods', numcells) );
                % overplot the index with transparency
                % first things first, edge detection:
                filt    = 'roberts';            % filter
                N_edges = edge( N, filt );
                % overplot the detected edges
                h       = imagesc( x_coords_all, obj.y_coords, repmat( ~N_edges, 1, 1, 3 ) );
                set( h, 'AlphaData', N_edges );
                set( gca, 'YDir', 'normal' );
            end

                
                
        end     % end function plotEz()
        
        
        function [obj] = plot_E_field_gui( obj )
            % plots all modes in a gui
            %
            %   Let's see what do i want
            %   Plot Ez, either with or without the phase and decay
            %   real imag amp or intensity
            %   with or without edges
            %   H field?
            %   with numcells
            %
            % Inputs:
            %   phi
            %       type: double, tensor
            %       desc: Field (envelope), dimensions y vs x vs mode #, where y = in
            %             plane transverse dimension and x = direction of propagation
            %   x
            %       type: double, vector
            %       desc: x coordinates
            %   y 
            %       type: double, vector
            %       desc: y coordinates
            %   k
            %       type: double, vector
            %       desc: prop. eigenvalues vs. mode #

            % Create a figure and axes
            f   = figure('Visible','off');
            ax  = axes('Units','pixels');

            % default settings
            mode_num            = 1;
            mode_comp           = 'real';       % 'real', 'imag', 'amplitude', 'intensity'
            phase_on            = true;         % determines whether the propagation phase is included
            index_overlay_on    = true;         % determines whether index overlay is on or off

            % number of modes
            n_modes = length( obj.k_vs_mode );

            % generate field
            numcells        = obj.numcells;
            nx              = round( numcells*obj.domain_size(2)/obj.dx );
            x_coords_all    = 0 : obj.dx : ( (nx-1) * obj.dx );
            % generate E field without phase (just the periodic envelope
            % function)
            Ez_no_phase     = repmat( obj.Phi_vs_mode, [ 1, numcells, 1 ] );
            % generate E field with phase (periodic envelope * exp( ikx ), including decay and growth)
            Ez_w_phase      = Ez_no_phase;
            % add phase in for loop for simplicity
            for ii = 1:n_modes

                % make phase
                phase_all = repmat( exp( 1i * obj.k_vs_mode(ii) * x_coords_all ), size(obj.Phi,1), 1 );

                % add in phase
                Ez_w_phase( :, :, ii ) = Ez_no_phase(:, :, ii) .* phase_all;

            end
        
            % repmat the index
            N_repmat = repmat( obj.N, 1, numcells );

            % initialize default field to be with phase. This is what gets
            % plotted
            Ez = Ez_w_phase;

            % plot the first mode by default
            plot_field( mode_num, mode_comp );

            % move axes over
            ax.Position(1) = ax.Position(1) + 100;

            % increase figure window size
            f.Position(4) = f.Position(4) + 25;
            f.Position(3) = f.Position(3) + 100;

            % create drop down menu labels
            nmodes          = length(obj.k_vs_mode);                % number of simulated modes
            labelstrings    = {};                                   % cell array of labels
            for ii = 1:nmodes
                labelstrings{end+1} = [ ' Mode ', num2str(ii) ];
            end

            % Create pop-up menu for selecting the mode
            popup_selectmode = uicontrol(  'Style', 'popup', ...
                                           'String', labelstrings, ...
                                           'Position', [20, ax.Position(4), 100, 50], ...
                                           'Callback', @selectmode );   

            % Create pop-up menu for selecting component to plot
            popup_selectcomponent = uicontrol( 'Style', 'popup', ...
                                               'String', {' Real', ' Imaginary', ' Amplitude', ' Intensity'}, ...
                                               'Position', [20, ax.Position(4) - 25, 100, 50], ...
                                               'Callback', @select_component );   

            % Create checkbox for selecting/deselecting phase
            checkbox_phase = uicontrol(    'Style', 'checkbox', ...
                                           'String', {'Phase on'}, ...
                                           'Value', phase_on, ...
                                           'Position', [20, ax.Position(4) - 50, 100, 30], ...
                                           'Callback', @toggle_phase );   
                                       
            % Create checkbox for turning grating index edges on and off
            checkbox_index_overlay = uicontrol(    'Style', 'checkbox', ...
                                                   'String', {'Index overlay on'}, ...
                                                   'Value', index_overlay_on, ...
                                                   'Position', [20, ax.Position(4) - 75, 100, 30], ...
                                                   'Callback', @toggle_index_overlay );   


            % Make figure visble after adding all components
            f.Visible = 'on';

            % function for selecting which mode to draw
            function selectmode( source, event )
                % selects which mode to draw

                mode_num     = source.Value;

                plot_field( mode_num, mode_comp );
            end

            % function for selecting which component to draw
            function select_component( source, event )
                % selects which mode to draw

                mode_comp    = source.String{ source.Value };
                mode_comp    = lower(mode_comp(2:end));                     % remove whitespace prefix 

                plot_field( mode_num, mode_comp );
            end

            % function for toggling phase on/off
            function toggle_phase( source, event )

                if source.Value == true
                    Ez = Ez_w_phase;
                else
                    Ez = Ez_no_phase;
                end 
                plot_field( mode_num, mode_comp );
            end
   
            % function for toggling index overlay on/off
            function toggle_index_overlay( source, event )

                if source.Value == true
                    index_overlay_on = true;
                else
                    index_overlay_on = false;
                end 
                plot_field( mode_num, mode_comp );
            end
            

            % function that plots the field
            function plot_field( mode_num, mode_comp )
                % depends on the variables:
                %   mode_num
                %   mode_comp
                % which are automatically updated by the other ui functions

                hold off;
                
                switch mode_comp

                    case 'real'
                        imagesc( x_coords_all, obj.y_coords, real(Ez(:,:,mode_num)) );
                        title( sprintf( 'Mode %i, real component', mode_num ));

                    case 'imaginary'
                        imagesc( x_coords_all, obj.y_coords, imag(Ez(:,:,mode_num)) );
                        title( sprintf( 'Mode %i, imaginary component', mode_num ));

                    case 'amplitude'
                        imagesc( x_coords_all, obj.y_coords, abs(Ez(:,:,mode_num)) );
                        title( sprintf( 'Mode %i, amplitude', mode_num ));

                    case 'intensity'
                        imagesc( x_coords_all, obj.y_coords, abs(Ez(:,:,mode_num)).^2 );
                        title( sprintf( 'Mode %i, intensity', mode_num ));

                end
                set( gca, 'ydir', 'normal' );
                colormap('redbluehilight');
                xlabel(['x (' obj.units.name ')']); ylabel(['y (' obj.units.name ')']);
                colorbar;

                if index_overlay_on == true
                    % overplot the index with transparency
                    
                    % first things first, edge detection:
                    filt    = 'roberts';            % filter
                    N_edges = edge( N_repmat, filt );
                    
                    % overplot the detected edges
                    hold on;
                    h       = imagesc( x_coords_all, obj.y_coords, repmat( ~N_edges, 1, 1, 3 ) );
                    set( h, 'AlphaData', N_edges );
                    set( gca, 'YDir', 'normal' );
                    
                end


            end     % end function plot_field()

        end     % end function plot_E_field_gui()
        
        
        function obj = runSimulation_v2_symm( obj, num_modes, BC, pml_options, guessk, OPTS )
            % Runs new, NEW HERMITIAN mode solver
            % exactly same as runSimulation but with new modesolver
            %
            % Description:
            %   Runs complex-k mode solver. Stores the mode with the most
            %   guided power. Calculates up/down power, directivity,
            %   angle of maximum radiation, and scattering strength.
            %
            % Inputs:
            %   num_modes
            %       type: integer
            %       desc: # of modes (max) to simulate
            %   BC
            %       type: integer
            %       desc: 0 for PEC, 1 for PMC
            %   pml_options
            %       type: array, double
            %       desc: 1x4 Array with the following elements:
            %               PML_options(1): PML in y direction (yes=1 or no=0)
            %               PML_options(2): length of PML layer in nm
            %               PML_options(3): strength of PML in the complex plane
            %               PML_options(4): PML polynomial order (1, 2, 3...)
            %   guessk
            %       type: scalar, double (can be complex)
            %       desc: guess k value. Works best when closest to desired
            %             mode. In units rad/'units'
            %   OPTS
            %       type: struct
            %       desc: optional options with the following fields
            %           'mode_to_overlap'
            %               type: matrix, double
            %               desc: mode to overlap
            %
            % Sets these properties:
            %   obj.k
            %       units rad/'units'
            %   obj.Phi
            %       dimensions y vs. x (x is dir. of propagation)
            %   obj.E_z
            %       field repeated, i can't remember why or if i use this
            %   obj.directivity
            %       up/down power ratio
            %   calls obj.calc_output_angle()
            %       this is a function, which probably figures out the
            %       output angle
            %   calls obj.calc_scattering_strength()
            %

            % default OPTS
            if nargin < 6
                OPTS = struct();
            end
            
            % spatial variables, in units nm
            nm          = 1e9;
            a           = obj.domain_size(2) * obj.units.scale * nm;
            lambda_nm   = obj.lambda * obj.units.scale * nm;
            dx_nm       = obj.dx * obj.units.scale * nm;
            guessk_nm   = guessk / ( obj.units.scale * nm );                % units rad/nm

            % store options
            obj.sim_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options, 'OPTS', OPTS );

            % set guessk if not entered
            if nargin < 5
                guessk = pi/(2*a);
            end
            
            % run solver
            k0_nm           = 2*pi/lambda_nm;
            [Phi_all, k_nm] = complexk_mode_solver_2D_PML_v2_symmetric( obj.N, ...
                                                       dx_nm, ...
                                                       k0_nm, ...
                                                       num_modes, ...
                                                       guessk_nm, ...
                                                       BC, ...
                                                       pml_options );
                                                   
            % re-scale k to units 'units'
            k_vs_mode = k_nm * nm * obj.units.scale;
            
            % save k and phi vs mode #
            obj.k_vs_mode   = k_vs_mode;
            obj.Phi_vs_mode = Phi_all;
                     
            % pick which mode to keep
            if isfield( OPTS, 'mode_to_overlap' )
                % pick mode with best overlap
                obj = obj.choose_mode( OPTS.mode_to_overlap );
            else
                % pick mode guided mode
                obj = obj.choose_mode(); 
            end
            
            % stitch together full e field, with the request number of
            % periods
            [obj, E_z]  = obj.stitch_E_field( obj.Phi, obj.k, obj.numcells );
            obj.E_z     = E_z;
            
            % for mode overlapping, stitch together single unit cell of E
            % field, using only real(k)
            [obj, E_z_for_overlap]  = obj.stitch_E_field( obj.Phi, real(obj.k), 1 );
            obj.E_z_for_overlap     = E_z_for_overlap;
            
            % pick slices of field to compute directivity, angle, etc.
            h_pml_d = round( pml_options(2)/obj.dy );                       % size of pml in discretizations
            y_up    = size( obj.E_z, 1 ) - h_pml_d - 1;
            y_down  = h_pml_d+2;
            
            % calculate output angle
            obj = obj.calc_output_angle( y_up, y_down );
            
            % calculate up/down directivity
            obj             = obj.calc_radiated_power();
            obj.directivity = obj.P_rad_up/obj.P_rad_down;
                    
            % calculate power scattering strength
            obj = obj.calc_scattering_strength();
            
            
        end     % end function runSimulation()
   
        
    end     % end methods
    
end     % end class


% -------------------------------------------------------------------------
% Legacy code
% -------------------------------------------------------------------------

%         function obj = runSimulation_old( obj, num_modes, BC, pml_options, guessk )
%             % LEGACY VERSION
%             % Runs old mode solver (jelena/mark)
%             %
%             % Description:
%             %   Runs complex-k mode solver. Stores the mode with the most
%             %   guided power. Calculates up/down power, directivity,
%             %   angle of maximum radiation, and scattering strength.
%             %
%             % Inputs:
%             %   num_modes
%             %       type: integer
%             %       desc: # of modes (max) to simulate
%             %   BC
%             %       type: integer
%             %       desc: 0 for PEC, 1 for PMC
%             %   pml_options
%             %       type: array, double
%             %       desc: 1x4 Array with the following elements:
%             %               PML_options(1): PML in y direction (yes=1 or no=0)
%             %               PML_options(2): length of PML layer in nm
%             %               PML_options(3): strength of PML in the complex plane
%             %               PML_options(4): PML polynomial order (1, 2, 3...)
% 
%             % spatial variables, in units nm
%             nm      = 1e9;
%             a       = obj.domain_size(2) * obj.units.scale * nm;
%             lambda  = obj.lambda * obj.units.scale * nm;
%             dx      = obj.dx * obj.units.scale * nm;
% 
%             % store options
%             obj.sim_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options );
% 
%             % set guessk if not entered
%             if nargin < 5
%                 guessk = pi/(2*a);
%             end
%             
%             % run solver
%             k0          = 2*pi/lambda;
%             [Phi_1D, k] = complexk_mode_solver_2D_PML_old( obj.N, ...
%                                                        dx, ...
%                                                        k0, ...
%                                                        num_modes, ...
%                                                        guessk, ...
%                                                        BC, ...
%                                                        pml_options );
%                                                    
%             % re-scale k
%             k = k * nm * obj.units.scale;
%             
%             % reshape field envelope ("Phi")
%             Phi_all = zeros( size(obj.N, 1), size(obj.N, 2), length(k) );  % stores all field envelope, dimensions y vs x vs k
%             for ii = 1:length(k)
%                 
%                 % reshape field
%                 Phi_k = reshape( Phi_1D(:,ii), fliplr( size(obj.N) ) ); % dimensions (x, y) where x = direction of propagation
%                 Phi_k = Phi_k.';
%                 
%                 % save in Phi_all
%                 Phi_all(:, :, ii) = Phi_k;  % dimensions (y, x, k)
%                 
%             end
%             
%             % DEBUG store temporary copies of k and phi_all b4 removing and
%             % sorting
%             temp_k_orig         = k;
%             temp_phi_all_orig   = Phi_all;
%             obj.debug.k_all     = temp_k_orig;
%             obj.debug.phi_all   = Phi_all;
% 
%             
%             % sort on guided power
%             guided_power  = zeros( size(k) );
%             total_power   = zeros( size(k) );
% 
%             % check to see if waveguide boundaries have been set yet
%             if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
%                 error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
%                         ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
%             end
%             
%             y_bot   = obj.wg_min_y;
%             y_top   = obj.wg_max_y;
%             y       = obj.y_coords;
%             
%             for ii = 1:length(k)
%                 % for each mode
% 
%                 % grab guided portion of field
%                 phi_guided = Phi_all( y >= y_bot & y <= y_top, :, ii );
%                 cur_phi    = Phi_all( :, :, ii );
% 
%                 % sum area of field
%                 guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
%                 total_power(ii)     = sum( abs( cur_phi(:) ).^2 );
% 
%             end
%             
%             % DEBUG storing the guided power
%             obj.debug.guided_power = guided_power;
%             
%             % keep mode with LEAST unguided power OR MOST guided power
%             [~, indx_k] = max( abs(guided_power./total_power) );        % most guided
%             k           = k(indx_k);
%             Phi         = Phi_all(:,:,indx_k);
%             
%             % check if k is backwards propagating
%             % if it is, then the field must be flipped.
%             if real(k) >=0 & imag(k) <= 0
%                 
%                 fprintf([ '\nMode found is backwards propagating (positive real k, negative imag k).\n', ...
%                           'Flipping the field and inverting the sign of k\n\n' ]);
%                 
%                 k   = -k;
%                 Phi = rot90(Phi, 2);    % equivalent to fliplr(flipud(Phi))
%                 
%             end
% 
%             % save wavevectors and field
%             obj.k   = k;
%             obj.Phi = Phi;
%             
%             % number of cells to repeat
%             numcells = obj.numcells;
% 
%             % stitch together e field, including the phase
%             x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
%             phase_all       = repmat( exp( 1i*k*x_coords_all ), size(Phi,1), 1 );
%             E_z             = repmat( Phi, 1, numcells ).*phase_all;
%             
%             % save E_z
%             obj.E_z         = E_z;
%             
%             % DEBUG plot the field
% %             obj.plotEz();
%             
%             % pick slices of field to compute directivity, angle, etc.
%             h_pml_d = round( pml_options(2)/obj.dy ); % size of pml in discretizations
%             y_up    = size(E_z,1) - h_pml_d - 1;
%             y_down  = h_pml_d+2;
%             
%             % calculate output angle
%             obj = obj.calc_output_angle( y_up, y_down );
%             
%             % calculate up/down directivity
%             obj             = obj.calc_radiated_power( y_up, y_down );
%             obj.directivity = obj.P_rad_up/obj.P_rad_down;
%                  
%             % calculate power scattering strength
%             obj = obj.calc_scattering_strength();
%             obj.debug.alpha_up_old    = imag(k) * obj.P_rad_up/( obj.P_rad_up + obj.P_rad_down );     % DEPRECATED upwards radiative loss
%             obj.debug.alpha_down_old  = imag(k) * obj.P_rad_down/( obj.P_rad_up + obj.P_rad_down);    % DEPRECATED downwards radiative loss
%             
%             
%         end     % end function runSimulation_old()
