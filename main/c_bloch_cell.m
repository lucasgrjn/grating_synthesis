classdef c_bloch_cell
    % Class for simulating a bloch cell, using FDFD complex-k solver
    %
    % SO i'm thinking that the current iteration of this is not that
    % useful, and probably overly complicated for grating design. But i may
    % come back to this for designing things other than gratings
    %
    % Authors: bohan zhang
    %
    % Description:
    %   - given input wavelength, solves for complex-k
    %   - currently only works for TE polarization
    %
    % outline of methods/properties
    %   - adds structures
    %   - runs complex-k fdfd solver
    %   - calculates things like directivity, angle, rad. eff
    %   - has dielectric profile
    %   - 
    
    properties
        
        OPTS;       % options
        N;          % index of refraction matrix
        constants;  % struct that holds various constants, i.e. speed of light, mu0, etc.
        p;          % TEMPORARY
        
        % struct that holds following fields
        %   'disc'   - discretization in both dimensions
        %   'units'  - saves units information (name = verbose units name,
        %               scale = scaling factor to convert to meters)
        %   'lambda' - wavelength, in units of 'units'
        %   'N'      - matrix of index of refraction
        domain;
        
    end
    
    methods
        
        % Constructor
        function obj = c_bloch_cell( varargin )
            % class constructor
            %
            % Inputs are name-value pairs
            % Inputs:
            % 
            %   'discretization'
            %       type: scalar, double
            %       desc: Discretization along both dimensions. Currently
            %             must be the same for both dimensions. Required
            %   
            %   'units'
            %       type: string
            %       desc: Spatial units, either 'm' for meters, 'mm' for
            %             millimetres, 'um' for micrometers, or 'nm' for
            %             nanometers.
            %             Defaults to 'um'
            %
            %   'lambda'
            %       type: double, scalar
            %       desc: simulation wavelength, in units of 'units'
            %
            %   'background_index'
            %       type: double, scalar
            %       desc: background index of simulation domain
            %             Defaults to 1.0 (air)
            %
            %   'domain_size'
            %       type: double, array
            %       desc: 1x2 array, [ y_size, z_size ] where y is
            %             transverse direction and z is propagation direction
            
            fprintf('\nConstructing new c_bloch_cell instance\n');
            
            % constant defs
            obj.constants.c     = 299792458;                                     % SI units [m/s]
            obj.constants.mu0   = 4e-7*pi;                                       % SI units [m kg/(s^2 A^2)]
            obj.constants.eps0  = 1/(obj.constants.c^2*obj.constants.mu0);       % SI units [s^4 A^2/(m^3 kg)]
            
            % list of required inputs, name-value pairs where value is
            % either the default or 'none' for no default (throws error if
            % not included in inputs)
            inputs = {  'discretization',   'none', ...
                        'units',            'um', ...
                        'lambda',           'none', ...
                        'background_index', 1.0, ...
                        'domain_size',      'none' ...
                     };     
                 
            % parse the inputs
            p = struct();
            for ii = 1:2:length(varargin)-1
                p.(varargin{ii}) = varargin{ii+1}; % = setfield( p, varargin{ii}, varargin{ii+1} );
            end
            
            % check existence of required inputs
            for ii = 1:2:( length(inputs)-1 )
               
                if ~isfield( p, inputs{ii} )
                    % this input was not found
                   
                    if ischar( inputs{ii+1} )
                        if strcmp( inputs{ii+1}, 'none' )
                            % this input has no default, throw error
                            error( 'Input ''%s'' was not set and requires a user-defined value. Try again.', inputs{ii} );
                        end
                    else
                        % this input has a default, set the default
                        fprintf( 'Input ''%s'' was not set, setting to default value ''%s''\n', inputs{ii}, inputs{ii+1} );
                        p.(inputs{ii}) = inputs{ii+1};
                    end
                    
                end
                
            end
            
            % -----------
            % Sort inputs
            obj.p           = p;          % DEBUG, temp code
            
            % domain inputs
            obj.domain.disc         = p.discretization;
            obj.domain.units.name   = p.units;
            switch( obj.domain.units.name )
                case 'm'
                    obj.domain.units.scale = 1;
                case 'mm'
                    obj.domain.units.scale = 1e-3;
                case 'um'
                    obj.domain.units.scale = 1e-6;
                case 'nm'
                    obj.domain.units.scale = 1e-9;
            end
            obj.domain.lambda   = p.lambda;
            obj.domain.size     = p.domain_size;
            % draw dielectric
            % TODO: check that integre # of discretization fits into domain
            % size
            % ALSO need to think about the grid locations....
            d               = obj.domain.disc;
            obj.domain.y    = 0:d:obj.domain.size(1)-d;
            obj.domain.z    = 0:d:obj.domain.size(2)-d;
            obj.domain.N    = p.background_index*ones( length(obj.domain.y), length(obj.domain.z) );
            
            fprintf('\n');
            
        end     % end constructor
        
        
        function obj = add_layer( obj ) %obj, min_x, height_x, indices, 
            % Adds horizontal layer(s) to the index
            layerFields = {'numLayers','centerZ','widthZ','minX','heightX',...
                'indices'};
            
        end
        
        
        function obj = run_complexk_solver( obj )
            % Not sure if this is useful, but wrapper function for running
            % the complex k mode solver
            
            
            
        end
        
        
    end
    
end

