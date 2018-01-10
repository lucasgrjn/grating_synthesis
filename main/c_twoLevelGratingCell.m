classdef c_twoLevelGratingCell
% authors: bohan zhang
% 
% encapsulates a 2 level grating cell simulation
% adapted from mark's code
%
% Notes to consider:
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
        k;              % complex k
        Phi;            % field envelope
        numcells;       % number of cells repeated in E_z
        E_z             % saves n_periods of the field, Phi(x,y)*exp(jkx) (z stands for z polarization)
        
        % radiation parameters
        directivity;        % ratio of power radiated up/power radiated down
        max_angle_up;       % Angle of maximum upwards output radiation
        max_angle_down;     % Angle of maximum downwards output radiation
        P_rad_down;         % saves power radiated down
        P_rad_up;           % saves power radiated up
        P_in;               % saves power entering the wg
        alpha_up;           % fractional radiative loss/unit length upwards
        alpha_down;         % fractional radiative loss/unit length downwards
        
        % DEBUG struct for holding temporary values that are useful during
        % debugging
        %   Current fields: k_all, phi_all, unguided_power, guided_power,
        %                   p_rad_up, Sx_up, Sx_down, Sy_up, Sy_down, Sx_in
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
            
            % create background dielectric
            obj.domain_size = p.domain_size;
            obj.x_coords    = 0 : obj.dx : obj.domain_size(2) - obj.dx;
            obj.y_coords    = 0 : obj.dy : obj.domain_size(1) - obj.dy;
            obj.N           = p.background_index * ones( length( obj.y_coords ), length( obj.x_coords ) );  % dimensions of y vs. x
            
            % check discretization fits in integer amt
            if abs(obj.x_coords(end) - (obj.domain_size(2) - obj.dx)) >= 1e-4
                % the discretization doesn't fit into x
                warning('dx doesn''t fit integer times into the x domain size.');
                
                % override the domain size
                fprintf('Overriding domain. Old domain x size was: %f. New domain x size is: %f\n\n',  obj.domain_size(2), obj.x_coords(end) + obj.dx);
                obj.domain_size(2) = obj.x_coords(end) + obj.dx;
            end
            if abs(obj.y_coords(end) - (obj.domain_size(1) - obj.dy)) >= 1e-4
                % the discretization doesn't fit into x
                warning('dy doesn''t fit integer times into the y domain size.');
                
                % override the domain size
                fprintf('Overriding domain. Old domain y size was: %f. New domain y size is: %f\n\n', obj.domain_size(1), obj.y_coords(end) + obj.dy);
                obj.domain_size(1) = obj.y_coords(end) + obj.dy;
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
            obj.N( y >= min_y & y < min_y + height_y, : ) = index; 
            
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
            obj.N( y >= min_y & y < min_y + height_y, x >= min_x & x < min_x + width_x ) = index; 
            
        end     % end function addRect()
        
        
        function obj = twoLevelBuilder( obj, wgs_min_y, wgs_thick, wgs_indx, ...
                                             wgs_duty_cycles, wgs_offsets )
            % Draws the two level waveguide grating
            % TODO: add option to include passivation layer (SiN)
            % sometime...
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
        
        
        function obj = runSimulation( obj, num_modes, BC, pml_options, guessk )
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

            % spatial variables, in units nm
            nm      = 1e9;
            a       = obj.domain_size(2) * obj.units.scale * nm;
            lambda  = obj.lambda * obj.units.scale * nm;
            dx      = obj.dx * obj.units.scale * nm;

            % store options
            obj.sim_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options );

            % set guessk if not entered
            if nargin < 5
                guessk = pi/(2*a);
            end
            
            % run solver
            k0          = 2*pi/lambda;
            [Phi_all, k] = complexk_mode_solver_2D_PML( obj.N, ...
                                                       dx, ...
                                                       k0, ...
                                                       num_modes, ...
                                                       guessk, ...
                                                       BC, ...
                                                       pml_options );
                                                   
            % re-scale k
            k = k * nm * obj.units.scale;
            
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
            
            % DEBUG store temporary copies of k and phi_all b4 removing and
            % sorting
            temp_k_orig         = k;
            temp_phi_all_orig   = Phi_all;
            obj.debug.k_all     = temp_k_orig;
            obj.debug.phi_all   = Phi_all;

            
            % sort on guided power
            guided_power  = zeros( size(k) );
            total_power   = zeros( size(k) );

            % check to see if waveguide boundaries have been set yet
            if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
                error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
                        ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
            end
            
            y_bot   = obj.wg_min_y;
            y_top   = obj.wg_max_y;
            y       = obj.y_coords;
            
            for ii = 1:length(k)
                % for each mode

                % grab guided portion of field
                phi_guided  = Phi_all( y >= y_bot & y <= y_top, :, ii );
                cur_phi     = Phi_all( :, :, ii );

                % sum area of field
                guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
                total_power(ii)     = sum( abs( cur_phi(:).^2 ) );

            end
            
            % DEBUG storing the guided power
            obj.debug.guided_power          = guided_power;
            obj.debug.guided_power_ratio    = guided_power./total_power;
            
            % keep mode with LEAST unguided power OR MOST guided power
            [~, indx_k] = max( abs(guided_power./total_power) );        % most guided
            k           = k(indx_k);
            Phi         = Phi_all(:,:,indx_k);
            
            % check if k is backwards propagating
            % if it is, then the field must be flipped.
            if real(k) > 0 && imag(k) < 0
                
                fprintf([ '\nMode found is backwards propagating (positive real k, negative imag k).\n', ...
                          'Re-running solver with flipped k\n\n' ]);

                % re-run solver
                [Phi, k] = complexk_mode_solver_2D_PML( obj.N, ...
                                                           dx, ...
                                                           k0, ...
                                                           1, ...
                                                           -k, ...
                                                           BC, ...
                                                           pml_options );

            end

            % save wavevectors and field
            obj.k   = k;
            obj.Phi = Phi;
            
            % number of cells to repeat
            numcells = obj.numcells;

            % stitch together e field, including the phase
            x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
            phase_all       = repmat( exp( 1i*k*x_coords_all ), size(Phi,1), 1 );
            E_z             = repmat( Phi, 1, numcells ).*phase_all;
            
            % save E_z
            obj.E_z         = E_z;
            
            % DEBUG plot the field
%             obj.plotEz();
            
            % pick slices of field to compute directivity, angle, etc.
            h_pml_d = round( pml_options(2)/obj.dy ); % size of pml in discretizations
            y_up    = size(E_z,1) - h_pml_d - 1;
            y_down  = h_pml_d+2;
            
            % calculate up/down directivity
            obj             = obj.calc_radiated_power( y_up, y_down );
            obj.directivity = obj.P_rad_up/obj.P_rad_down;
            
            % calculate output angle
            obj = obj.calc_output_angle( y_up, y_down );
            
            % calculate power scattering strength
            obj                       = obj.calc_scattering_strength();
            obj.debug.alpha_up_old    = imag(k) * obj.P_rad_up/( obj.P_rad_up + obj.P_rad_down );     % DEPRECATED upwards radiative loss
            obj.debug.alpha_down_old  = imag(k) * obj.P_rad_down/( obj.P_rad_up + obj.P_rad_down);    % DEPRECATED downwards radiative loss
            
            
        end     % end function runSimulation()
        
        
        function obj = runSimulation_old( obj, num_modes, BC, pml_options, guessk )
            % LEGACY VERSION
            % Runs old mode solver (jelena/mark)
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

            % spatial variables, in units nm
            nm      = 1e9;
            a       = obj.domain_size(2) * obj.units.scale * nm;
            lambda  = obj.lambda * obj.units.scale * nm;
            dx      = obj.dx * obj.units.scale * nm;

            % store options
            obj.sim_opts = struct( 'num_modes', num_modes, 'BC', BC, 'pml_options', pml_options );

            % set guessk if not entered
            if nargin < 5
                guessk = pi/(2*a);
            end
            
            % run solver
            k0          = 2*pi/lambda;
            [Phi_1D, k] = complexk_mode_solver_2D_PML_old( obj.N, ...
                                                       dx, ...
                                                       k0, ...
                                                       num_modes, ...
                                                       guessk, ...
                                                       BC, ...
                                                       pml_options );
                                                   
            % re-scale k
            k = k * nm * obj.units.scale;
            
            % reshape field envelope ("Phi")
            Phi_all = zeros( size(obj.N, 1), size(obj.N, 2), length(k) );  % stores all field envelope, dimensions y vs x vs k
            for ii = 1:length(k)
                
                % reshape field
                Phi_k = reshape( Phi_1D(:,ii), fliplr( size(obj.N) ) ); % dimensions (x, y) where x = direction of propagation
                Phi_k = Phi_k.';
                
                % save in Phi_all
                Phi_all(:, :, ii) = Phi_k;  % dimensions (y, x, k)
                
            end
            
            % DEBUG store temporary copies of k and phi_all b4 removing and
            % sorting
            temp_k_orig         = k;
            temp_phi_all_orig   = Phi_all;
            obj.debug.k_all     = temp_k_orig;
            obj.debug.phi_all   = Phi_all;

            
            % sort on guided power
            guided_power  = zeros( size(k) );
            total_power   = zeros( size(k) );

            % check to see if waveguide boundaries have been set yet
            if isempty(obj.wg_min_y) || isempty(obj.wg_max_y)
                error(['Waveguide boundaries have not been set yet. You must set the waveguide boundaries by either calling' ...
                        ' "twoLevelBuilder()" or by setting the "wg_min_y" AND "wg_max_y" object properties yourself.']);
            end
            
            y_bot   = obj.wg_min_y;
            y_top   = obj.wg_max_y;
            y       = obj.y_coords;
            
            for ii = 1:length(k)
                % for each mode

                % grab guided portion of field
                phi_guided = Phi_all( y >= y_bot & y <= y_top, :, ii );
                cur_phi    = Phi_all( :, :, ii );

                % sum area of field
                guided_power(ii)    = sum( abs( phi_guided(:) ).^2 );
                total_power(ii)     = sum( abs( cur_phi(:) ).^2 );

            end
            
            % DEBUG storing the guided power
            obj.debug.guided_power = guided_power;
            
            % keep mode with LEAST unguided power OR MOST guided power
            [~, indx_k] = max( abs(guided_power./total_power) );        % most guided
            k           = k(indx_k);
            Phi         = Phi_all(:,:,indx_k);
            
            % check if k is backwards propagating
            % if it is, then the field must be flipped.
            if real(k) >=0 & imag(k) <= 0
                
                fprintf([ '\nMode found is backwards propagating (positive real k, negative imag k).\n', ...
                          'Flipping the field and inverting the sign of k\n\n' ]);
                
                k   = -k;
                Phi = rot90(Phi, 2);    % equivalent to fliplr(flipud(Phi))
                
            end

            % save wavevectors and field
            obj.k   = k;
            obj.Phi = Phi;
            
            % number of cells to repeat
            numcells = obj.numcells;

            % stitch together e field, including the phase
            x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
            phase_all       = repmat( exp( 1i*k*x_coords_all ), size(Phi,1), 1 );
            E_z             = repmat( Phi, 1, numcells ).*phase_all;
            
            % save E_z
            obj.E_z         = E_z;
            
            % DEBUG plot the field
%             obj.plotEz();
            
            % pick slices of field to compute directivity, angle, etc.
            h_pml_d = round( pml_options(2)/obj.dy ); % size of pml in discretizations
            y_up    = size(E_z,1) - h_pml_d - 1;
            y_down  = h_pml_d+2;
            
            % calculate up/down directivity
            obj             = obj.calc_radiated_power( y_up, y_down );
            obj.directivity = obj.P_rad_up/obj.P_rad_down;
            
            % calculate output angle
            obj = obj.calc_output_angle( y_up, y_down );
            
            % calculate power scattering strength
            obj = obj.calc_scattering_strength();
            obj.debug.alpha_up_old    = imag(k) * obj.P_rad_up/( obj.P_rad_up + obj.P_rad_down );     % DEPRECATED upwards radiative loss
            obj.debug.alpha_down_old  = imag(k) * obj.P_rad_down/( obj.P_rad_up + obj.P_rad_down);    % DEPRECATED downwards radiative loss
            
            
        end     % end function runSimulation_old()
        
        
        function obj = calc_radiated_power(obj, y_up, y_down)
            % Calculates the radiated power in both the upwards and
            % downwards directions
            % Only run AFTER mode solver has been run, since this function
            % depends on E_z and k
            %
            % Inputs:
            %   y_up    - location of top slice of E field to take
            %   y_down  - location of bot slice of E field to take
            %
            % Uses these object properties:
            %   E_z
            %   k
            %
            % Sets these object properties:
            %   P_rad_down      - saves power radiated down
            %   P_rad_up        - saves power radiated up
            
            % define constants
            mu0     = 4*pi*1e-7;                % units of H/m
            mu0     = mu0 * obj.units.scale;    % units of H/(units)
            c       = 3e8;                      % units of m/s
            c       = c/obj.units.scale;        % units of (units)/s
            omega0 	= 2*pi*c/obj.lambda;     % units of rad*(Units/s)/units = rad/s
            
            % grab k and E_z
            E_z = obj.E_z;
            k   = obj.k; 
            
            % Calculate power radiated up
            
            % calculate H_x
            H_x_up  = 1/(1i*omega0*mu0) .* ( E_z(y_up-1,:) - E_z(y_up+1,:) )/(2*obj.dy);
            
            % calculate H_y
            % calculate full H_y, see notes for derivation
            H_y_up1 = (1i/(omega0*mu0)) * ( E_z(y_up, 3:2:end ) - E_z(y_up, 1:2:end-2 ) )/(2*obj.dx);   % dx, term staggered 1
            H_y_up2 = (1i/(omega0*mu0)) * ( E_z(y_up, 4:2:end ) - E_z(y_up, 2:2:end-2 ) )/(2*obj.dx);   % dx, term staggered 2
            % interleave two arrays
            H_y_up                                  = zeros( size( [H_y_up1, H_y_up2] ) );
            H_y_up( ( 1:length(H_y_up1) )*2 - 1)    = H_y_up1;
            H_y_up( ( 1:length(H_y_up2) )*2 )       = H_y_up2;

            % power radiated up
            Sy_up = real( H_x_up(:)' .* E_z(y_up,:) );                          % Sy = real( Ez Hx* )
%             Sx_up = real( (-1) * H_y_up(:)' .* E_z( y_up, 2:end-1 ) );          % Sx = real( -Ez Hy* )
%             P_rad_up_Hx_only    = sum( abs(Sy_up) )*obj.dx;                                    % only using Sy component
%             P_rad_up      = sum( sqrt( Sy_up( 2:end-1 ).^2 + Sx_up.^2 ) )*obj.dx;        % using both Sx and Sy compmonents
            P_rad_up      = sum( abs(Sy_up) )*obj.dx;        % using both Sx and Sy compmonents

            % DEBUG save p rad up
%             obj.debug.P_rad_up_Hx_only  = P_rad_up_Hx_only;
%             obj.debug.P_rad_up_Hx_Hy    = P_rad_up_Hx_Hy;
            
            % calculate power radiated down
            
            % calculate H_x
            H_x_down  = 1/(1i*omega0*mu0) .* ( E_z(y_down-1,:) - E_z(y_down+1,:) )/(2*obj.dy);
            
            % calculate H_y
            % calculate full H_y, see notes for derivation
            H_y_down1 = (1i/(omega0*mu0)) * ( E_z(y_down, 3:2:end ) - E_z(y_down, 1:2:end-2 ) )/(2*obj.dx);   % dx, term staggered 1
            H_y_down2 = (1i/(omega0*mu0)) * ( E_z(y_down, 4:2:end ) - E_z(y_down, 2:2:end-2 ) )/(2*obj.dx);   % dx, term staggered 2
            % interleave two arrays
            H_y_down                                    = zeros( size( [H_y_down1, H_y_down2] ) );
            H_y_down( ( 1:length(H_y_down1) )*2 - 1)    = H_y_down1;
            H_y_down( ( 1:length(H_y_down2) )*2 )       = H_y_down2;

            % power radiated down
            Sy_down = real( H_x_down(:)' .* E_z(y_down,:) );                % Sy = real( Ez Hx* )
%             Sx_down = real( -H_y_down(:)' .* E_z(y_down, 2:end-1 ) );       % Sx = real( -Ez Hy* )
%             P_rad_down_Hx_only  = sum( abs(Sy_down) );                      % only using Sy component
            P_rad_down    = sum( abs(Sy_down) )*obj.dx;     % using both Sx and Sy compmonents

            % DEBUG save p rad down
%             obj.debug.P_rad_down_Hx_only  = P_rad_down_Hx_only;
%             obj.debug.P_rad_down_Hx_Hy    = P_rad_down_Hx_Hy;

            % save power radiated up and down
%             P_rad_up   = P_rad_up_Hx_Hy;
%             P_rad_down = P_rad_down_Hx_Hy;
            obj.P_rad_down    = P_rad_down;
            obj.P_rad_up      = P_rad_up;
            
            % DEBUG save Sx and Sy
%             obj.debug.Sx_up     = Sx_up;
%             obj.debug.Sx_down   = Sx_down;
            obj.debug.Sy_up     = Sy_up;
            obj.debug.Sy_down   = Sy_down;
            
        end     % end function calc_radiated_power()
        
        
        function obj = calc_output_angle(obj, y_up, y_down)
            % Calculates output angle of radiation
            % To be run only AFTER mode solver has been run
            %
            % Inputs:
            %   y_up    - location of top slice of E field to take
            %   y_down  - location of bot slice of E field to take
            
            % load stuff
            E_z = obj.E_z;
            f0  = 1/obj.lambda;     % spatial freq
            
            % grab slices
            E_z_up          = E_z( y_up, : );
            E_z_down        = E_z( y_down, : );
            
            % do some zero padding
            pad_factor  = 10;
            E_z_up      = [ zeros( 1, pad_factor*length(E_z_up) ), E_z_up, zeros( 1, pad_factor*length(E_z_up) ) ];
            E_z_down    = [ zeros( 1, pad_factor*length(E_z_down) ), E_z_down, zeros( 1, pad_factor*length(E_z_down) ) ];
            
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
            
            % NEW: H y in , on half step, using entire E_z slice
            H_y_in = (1i/(omega0*mu0)) * ( E_z( :, 2 ) - E_z( :, 1 ) )/(obj.dx);   % dy, term staggered 1
            
            % power in (at left edge)
%             Sy_in = real( conj(H_x_in(:)) .* E_z( 2:end-1, 2 ) );           % Sy = real( Ez Hx* )
            Sx_in = real( -1 * conj(H_y_in(:)) .* E_z( :, 2 ) );      % Sx = real( -Ez Hy* )
%             P_in  = sum( sqrt( Sy_in.^2 + Sx_in.^2 ) )*obj.dy;              % using both Sx and Sy compmonents
%             % NEW use just Sx
            P_in = sum( Sx_in(:) )*obj.dy;

%             % DEBUG let H_y_in be on half step'
%             H_y_in_half = (1i/(omega0*mu0)) * ( E_z( 2:end-1, 2 ) - E_z( 2:end-1, 1 ) )/(obj.dx);   % dy, term staggered 1
%             Sx_in_athalf = real( -1 * conj(H_y_in_half(:)) .* E_z( 2:end-1, 1 ) );                         % Sx = real( -Ez Hy* )
%             P_in_athalf  = sum( sqrt( Sy_in.^2 + Sx_in_athalf.^2 ) )*obj.dy;              % using both Sx and Sy compmonents
            
            % DEBUG save Sx_in
            obj.debug.Sx_in = Sx_in;

            % save input power
            obj.P_in = P_in;
            % TEMP DEBUG
%             obj.P_in = P_in_athalf;
            
            % calculate alpha efficiency (in 1/units)
            period          = obj.domain_size(2);
            numcells        = obj.numcells;
%             obj.alpha_up    = ( -1/( 2*numcells*period ) ) * log( obj.P_rad_up/P_in );
%             obj.alpha_down  = ( -1/( 2*numcells*period ) ) * log( obj.P_rad_down/P_in );
            obj.alpha_up    = ( -1/( 2*numcells*period ) ) * log( 1 - obj.P_rad_up/P_in );
            obj.alpha_down  = ( -1/( 2*numcells*period ) ) * log( 1 - obj.P_rad_down/P_in );
%             obj.alpha_tot   = obj.alpha_up + obj.alpha_down;
                     
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
        
        
        function plotEz_w_edges(obj)
            % plots electric field, this could probably be part of a gui
            % with grating geometry edges overlaid
            
            numcells = obj.numcells;
            E_z      = obj.E_z;
            
            % repmat the index
            N = repmat( obj.N, 1, numcells );
            
            % define x coordinates
            x_coords_all    = 0 : obj.dx : numcells*obj.domain_size(2)-obj.dx;
            
                        % plot the field, real
            figure;
            imagesc( x_coords_all, obj.y_coords, real(E_z) ); hold on;
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

            % plot the field, abs
            figure;
            imagesc( x_coords_all, obj.y_coords, abs(E_z) ); hold on;
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
            
        end     % end function plotEz()
        
    end
    
end

