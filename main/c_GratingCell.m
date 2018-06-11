classdef c_GratingCell < c_bloch_cell
% Class for simulating a grating cell, using FDFD complex-k solver
%
% Authors: bohan zhang
%
% Description:
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
% Example usage:
    
    properties
        
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
        
    end     % end properties
    
    
    methods
        
        % -----------------------------------------------------------------
        % Constructor
        
        function obj = c_GratingCell( varargin )
            
            % call bloch cell constructor
            obj = obj@c_bloch_cell(varargin);
            
        end     % end constructor
        
        % -----------------------------------------------------------------
        % Running simulations
        
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
            
            % first run the bloch cell simulation function
            obj = runSimulation@c_bloch_cell( num_modes, BC, pml_options, guessk, OPTS );
            
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
            % i don't even remember what this is for
            % looks like it's used in calc_scattering_strength
            [obj, obj.E_z]  = obj.stitch_E_field( obj.Phi, obj.k, obj.numcells );
            
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
            
        end     % end runSimulation
        
        
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
        
    end     % end methods
    
    
end

