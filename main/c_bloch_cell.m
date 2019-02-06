%> Class for simulating a bloch cell, using FDFD complex-k solver
%>
%> Authors: bohan zhang
%>
%> Description:
%>
%> Notes:
%>   - propagation direction coordinate is 'x'. 
%>   - transverse (in plane) direction coordinate is 'y'.
%>   - transverse (out of plane) direction is 'z'.
%>
%> Inputs to constructor: \n 
%>   Inputs are name-value pairs: \n
%>   'discretization' \n
%>       type: double, scalar or 1x2 vector \n
%>       desc: discretization along x and y, in units of 'units' \n
%>             if scalar, then dx = dy = discretization \n
%>             if vector, then discretization = [ dy dx ] \n
%>
%>   'units' \n
%>       type: string \n
%>       desc: name and scaling of spatial units, supports 'm' \n
%>             (meters), 'mm' (millimeters), 'um' (microns), 'nm' \n
%>             (nanometers) \n
%>
%>   'lambda' \n
%>       type: double, scalar \n
%>       desc: wavelength to solve at, in units 'units' \n
%>
%>   'background_index' \n
%>       type: double, scalar \n
%>       desc: value of background index \n
%>
%>   'domain_size' \n
%>       type: 1x2 array, double \n
%>       desc: domain size, [ y height, x length ] \n
%>
%>   'num_cells' \n
%>       type: integer, scalar \n
%>       desc: OPTIONAL: pick the number of grating cells to repeat when \n
%>             drawing/calculating Ez \n
%>
%> Example usage:
classdef c_bloch_cell
    
    properties
        
        %> index profile
        N;              
        %> discretization in x (dir of propagation)
        dx;                 
        %> discretization in y (transverse direction)
        dy;
        %> name and scaling of spatial units, supports 'm', 'mm', 'um', 'nm'
        units;
        %> free space wavelength, units 'units'
        lambda;
        %> [ max_y, max_x ], where y = transverse and x = direction of propagation
        domain_size;
        %> vector of x coordinates (dir of prop)
        x_coords;
        %> vector of y coords (transverse dir)
        y_coords;           
        
        %> struct that holds simulation options \n
        %> current options are: \n
        %>  'num_modes'     - number of modes \n
        %>   'BC'            - boundary conditions, 0 for pec 1 for pmc i think \n
        %>   'pml_options'   - pml options, see complexk solver for details \n
        %>                     PML_options(1): PML in y direction (yes=1 or no=0) \n
        %>                     PML_options(2): length of PML layer in nm \n
        %>                     PML_options(3): strength of PML in the complex plane \n
        %>                     PML_options(4): PML polynomial order (1, 2, 3...) \n
        sim_opts;
        
        %> mode characteristics
        %> complex k, units rad/'units'
        k;                  
        %> field envelope
        Phi;                
        %> number of cells repeated in E_z
        numcells;           
        %> saves n_periods of the field, Phi(x,y)*exp(jkx) (z stands for z polarization)
        E_z                 
        %> k vs mode #
        k_vs_mode;          
        %> Phi vs mode #, dimensions y vs. x vs mode #
        Phi_vs_mode;        
        %> which mode was chosen
        chosen_mode_num;    
        %> Field for mode overlapping. One unit cell with only real(k) in the phase
        E_z_for_overlap;    
        
        %> DEBUG struct for holding temporary values that are useful during
        %> debugging
        %>   Current fields: k_all, phi_all, unguided_power, guided_power,
        %>                   p_rad_up, Sx_up, Sx_down, Sy_up, Sy_down, Sx_in
        %>                   P_rad_up_onecell, P_rad_down_onecell,
        %>                   Sx, Sy, P_per_y_slice, P_per_x_slice
        %>   Some of the above fields may have been removed.
        debug;
        
    end     % end properties
    
    
    methods
        
        % -----------------------------------------------------------------
        %> Constructor
        function obj = c_bloch_cell( varargin )
            
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
            
            % new version of creating background dielectric that rounds
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
            
            % set number of cells, not sure if used anymore
            obj.numcells = p.numcells;
            
        end     % end constructor
        
        % -----------------------------------------------------------------
        % Drawing functions
        
        % -----------------
        %> Draws a horizontal layer of dielectric
        %>
        %> Inputs:
        %>   min_y
        %>       Desc: scalar double, minimum y
        %>   height_y
        %>       Desc: scalar double, height/thickness of layer
        %>   index
        %>       Desc: scalar double, index of refraction of layer
        function obj = addLayer( obj, min_y, height_y, index )

            
            y = obj.y_coords;
            
            % fill in the layer
            % cut off before end of grid
            obj.N( y > (min_y - obj.dy/2) & y < (min_y + height_y - obj.dy/2), : ) = index; 
            
        end     % end function addLayer()
        
        % -----------------
        %> Draws a rectangle
        %>
        %> Inputs:
        %>   min_x
        %>       Desc: scalar double, left edge of rectangle
        %>   min_y
        %>       Desc: scalar double, bottom edge of rectangle
        %>   width_x
        %>       Desc: scalar double, width of rectangle
        %>   height_y
        %>       Desc: scalar double, height of rectangle
        %>   index
        %>       Desc: index of refraction
        function obj = addRect( obj, min_x, min_y, width_x, height_y, index )
            
            
            x = obj.x_coords;
            y = obj.y_coords;
            
            % fill in the rect
            % cut off before end of grid
            obj.N( y > (min_y - obj.dy/2) & y < (min_y + height_y - obj.dy/2), ...
                    x > (min_x - obj.dx/2) & x < (min_x + width_x - obj.dx/2) ) = index; 
%             % cut off after end of grid
%             obj.N( y > (min_y - obj.dy/2) & y < (min_y + height_y - obj.dy/2), ...
%                     x > (min_x - obj.dx/2) & x < (min_x + width_x - obj.dx/2) ) = index; 
%             % cut off right before grid
%             obj.N( y > (min_y - obj.dy/2) & y < (min_y + height_y), ...
%                     x > (min_x - obj.dx/2) & x < (min_x + width_x) ) = index; 
%             % cut off on grid
%             obj.N( y > (min_y - obj.dy/2) & y <= (min_y + height_y), ...
%                     x > (min_x - obj.dx/2) & x <= (min_x + width_x) ) = index; 
%             % cut off on grid, both ends
%             obj.N( y >= (min_y) & y <= (min_y + height_y), ...
%                     x >= (min_x) & x <= (min_x + width_x) ) = index; 
%             % cut off on grid on left, right before grid on right
%             obj.N( y >= (min_y) & y < (min_y + height_y), ...
%                     x >= (min_x) & x < (min_x + width_x) ) = index; 
            
        end     % end function addRect()
        
        % -----------------------------------------------------------------
        % Running simulations
        
        % -----------------
        %> Runs new mode solver
        %>
        %> Description:
        %>   Runs complex-k mode solver. Stores the mode with the most
        %>   guided power. Calculates up/down power, directivity,
        %>   angle of maximum radiation, and scattering strength.
        %>
        %> Inputs:
        %>   num_modes
        %>       type: integer
        %>       desc: # of modes (max) to simulate
        %>   BC
        %>       type: integer
        %>       desc: 0 for PEC, 1 for PMC
        %>   pml_options
        %>       type: array, double
        %>       desc: 1x4 Array with the following elements:
        %>               PML_options(1): PML in y direction (yes=1 or no=0)
        %>               PML_options(2): length of PML layer in nm
        %>               PML_options(3): strength of PML in the complex plane
        %>               PML_options(4): PML polynomial order (1, 2, 3...)
        %>   guessk
        %>       type: scalar, double (can be complex)
        %>       desc: guess k value. Works best when closest to desired
        %>             mode. In units rad/'units'
        %>   OPTS
        %>       type: struct
        %>       desc: optional options with the following fields
        %>           'mode_to_overlap'
        %>               type: matrix, double
        %>               desc: mode to overlap
        %>
        %> Sets these properties:
        %>   obj.k
        %>       units rad/'units'
        %>   obj.Phi
        %>       dimensions y vs. x (x is dir. of propagation)
        %>   obj.E_z
        %>       field repeated, i can't remember why or if i use this
        %>   obj.directivity
        %>       up/down power ratio
        function obj = runSimulation( obj, num_modes, BC, pml_options, guessk, OPTS )

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
            
            % The stuff below will be moved to grating cell
%             % pick which mode to keep
%             if isfield( OPTS, 'mode_to_overlap' )
%                 % pick mode with best overlap
%                 obj = obj.choose_mode( OPTS.mode_to_overlap );
%             else
%                 % pick mode guided mode
%                 obj = obj.choose_mode(); 
%             end
            
%             % stitch together full e field, with the request number of
%             % periods
%             [obj, E_z]  = obj.stitch_E_field( obj.Phi, obj.k, obj.numcells );
%             obj.E_z     = E_z;
            
%             % for mode overlapping, stitch together single unit cell of E
%             % field, using only real(k)
%             [obj, E_z_for_overlap]  = obj.stitch_E_field( obj.Phi, real(obj.k), 1 );
%             obj.E_z_for_overlap     = E_z_for_overlap;
            
%             % pick slices of field to compute directivity, angle, etc.
%             h_pml_d = round( pml_options(2)/obj.dy );                       % size of pml in discretizations
%             y_up    = size( obj.E_z, 1 ) - h_pml_d - 1;
%             y_down  = h_pml_d+2;
            
%             % calculate output angle
%             obj = obj.calc_output_angle( y_up, y_down );
%             
%             % calculate up/down directivity
%             obj             = obj.calc_radiated_power();
%             obj.directivity = obj.P_rad_up/obj.P_rad_down;
%                     
%             % calculate power scattering strength
%             obj = obj.calc_scattering_strength();
                
        end     % end function runSimulation()
        
        
        % -----------------
        %> Function that chooses which mode becomes the accepted mode
        %>
        %> Inputs:
        %>   mode_to_overlap
        %>       type: matrix, double
        %>       desc: This function will choose the mode with the closest overlap.
        %>             The mode_to_overlap doesn't have to have the same
        %>             size as the domain
        %>
        %> Sets these properties:
        %>   obj.k
        %>       prop. k of chosen mode
        %>   obj.Phi
        %>       field envelope of chosen mode
        %>   obj.chosen_mode_num
        %>       index of chosen mode (chosen mode k = obj.k_vs_mode(
        %>       obj.chosen_mode_num))
        function obj = choose_mode( obj, mode_to_overlap )

            % run overlaps
            [obj, max_overlaps] = obj.calc_mode_overlaps( mode_to_overlap );

            % keep mode with highest overlap
            [~, indx_k]         = max( max_overlaps );        
            obj.k               = obj.k_vs_mode(indx_k);
            obj.Phi             = obj.Phi_vs_mode(:,:,indx_k);
            obj.chosen_mode_num = indx_k;
            
        end     % end function choose_mode()
        
        
        % -----------------
        %> Takes in a mode and overlaps it with the modes that were
        %> solved in this object
        %> Overlap is performed using a cross correlation
        %>
        %> Inputs:
        %>   mode_to_overlap
        %>       type: matrix, double
        %>       desc: field to overlap with. Remove the exponential
        %>             term from it
        %>
        %> Outputs:
        %>   max_overlaps
        %>       type: vector, double
        %>       desc: Max mode overlap vs. mode #   
        function [obj, max_overlaps] = calc_mode_overlaps( obj, mode_to_overlap )  
            
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
        
        
        % -----------------
        %> Stitches the E field together from the phase and envelope
        %>
        %> Inputs:
        %>   Phi
        %>       type: matrix, double
        %>       desc: field envelope
        %>   k
        %>       type: scalar, double
        %>       desc: propagation constant
        %>   num_cells
        %>       type: scalar, int
        %>       desc: number of cells to repeat
        function [obj, E_z] = stitch_E_field( obj, Phi, k, num_cells )

            
            % stitch together e field, including the phase
            nx              = round( num_cells*obj.domain_size(2)/obj.dx );
            x_coords_all    = 0 : obj.dx : ( (nx-1) * obj.dx );
            phase_all       = repmat( exp( 1i*k*x_coords_all ), size(Phi,1), 1 );
            E_z             = repmat( Phi, 1, num_cells ).*phase_all;
            
        end     % end function stitch_E_field()
        
        % -----------------
        %> Circularly shifts the index by the shift length
        %> Shift length can be positive or negative
        %> Ideally the shift length is an integer multiple of obj.dx
        %>
        %> Inputs:
        %>   shift_length_x
        %>       type: double, scalar
        %>       desc: length to shift by, can be either positive or
        %>             negative value
        function obj = shift_index_circ( obj, shift_length_x )
            
            % calc number of units to shift by
            nx      = round( shift_length_x/obj.dx );
            obj.N   = circshift( obj.N, nx, 2 );
            
        end     % end function shift_index_circ()
        
        
        % -----------------
        %> Plots the index distribution
        function plotIndex(obj)
                        
            figure;
            imagesc( obj.x_coords, obj.y_coords, obj.N );
            colorbar;
            set( gca, 'YDir', 'normal' );
            xlabel(['x (', obj.units.name, ')']);
            ylabel(['y (', obj.units.name, ')']);
            title('Plot of index distribution');
            
        end     % end function plotIndex()
        
        
        % -----------------
        %> plots all modes in a gui
        %>
        %>   Plot Ez, either with or without the phase and decay
        %>   real imag amp or intensity
        %>   with or without edges
        %>   H field?
        %>   with numcells
        %>
        %> Inputs:
        %>   phi
        %>       type: double, tensor
        %>       desc: Field (envelope), dimensions y vs x vs mode #, where y = in
        %>             plane transverse dimension and x = direction of propagation
        %>   x
        %>       type: double, vector
        %>       desc: x coordinates
        %>   y 
        %>       type: double, vector
        %>       desc: y coordinates
        %>   k
        %>       type: double, vector
        %>       desc: prop. eigenvalues vs. mode #
        function [obj] = plot_E_field_gui( obj )

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

            % Create textbox for displaying the mode's eigenvalue
            textbox_eigenval = uicontrol( 'Style', 'text', ...
                                          'String', 'hello', ...
                                          'Position', [20, ax.Position(4) - 100, 100, 30] );
                                           
                                               

            % Make figure visble after adding all components
            f.Visible = 'on';

            %> function for selecting which mode to draw
            function selectmode( source, event )
                % selects which mode to draw

                mode_num     = source.Value;

                plot_field( mode_num, mode_comp );
            end

            %> function for selecting which component to draw
            function select_component( source, event )
                % selects which mode to draw

                mode_comp    = source.String{ source.Value };
                mode_comp    = lower(mode_comp(2:end));                     % remove whitespace prefix 

                plot_field( mode_num, mode_comp );
            end

            %> function for toggling phase on/off
            function toggle_phase( source, event )

                if source.Value == true
                    Ez = Ez_w_phase;
                else
                    Ez = Ez_no_phase;
                end 
                plot_field( mode_num, mode_comp );
            end
   
            %> function for toggling index overlay on/off
            function toggle_index_overlay( source, event )

                if source.Value == true
                    index_overlay_on = true;
                else
                    index_overlay_on = false;
                end 
                plot_field( mode_num, mode_comp );
            end
            

            %> function that plots the field
            %> depends on the variables:
            %>   mode_num
            %>   mode_comp
            %> which are automatically updated by the other ui functions
            function plot_field( mode_num, mode_comp )

                hold off;
                
                % E field to plot
                Ez_to_plot = Ez(:,:,mode_num);
                Ez_to_plot = Ez_to_plot./max(abs(Ez_to_plot(:)));
                
                switch mode_comp

                    case 'real'
                        imagesc( x_coords_all, obj.y_coords, real(Ez_to_plot) );
                        title( sprintf( 'Mode %i, real component', mode_num ));

                    case 'imaginary'
                        imagesc( x_coords_all, obj.y_coords, imag(Ez_to_plot) );
                        title( sprintf( 'Mode %i, imaginary component', mode_num ));

                    case 'amplitude'
                        imagesc( x_coords_all, obj.y_coords, abs(Ez_to_plot) );
                        title( sprintf( 'Mode %i, amplitude', mode_num ));

                    case 'intensity'
                        imagesc( x_coords_all, obj.y_coords, abs(Ez_to_plot).^2 );
                        title( sprintf( 'Mode %i, intensity', mode_num ));

                end
                set( gca, 'ydir', 'normal' );
                colormap('redbluehilight');
                caxis( [-1,1] .* max(abs(Ez_to_plot(:))) );                 % get the color limits, center white = 0
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
        
    end     % end methods
    
    
end     % end classdef

