classdef c_layer
    % Layer primitive
    % Stores y location, thickness, name, material
    
    properties
        
        % minimum y position
        min_y; 
        
        % layer thickness
        height_y;
        
        % index of refraction
        index;
        
        % name of this layer
        name;
        
    end
    
    methods
        
        function obj = c_layer( min_y, height_y, index, name )
            % Constructor
            %
            % Inputs:
            %   min_y
            %   height_y
            %   index
            %   name
            
            obj.min_y       = min_y;
            obj.height_y    = height_y;
            obj.index       = index;
            obj.name        = name;
            
            
        end     % end constructor
        
        
        function new_diel = draw( obj, diel )
            
        end
        
    end
    
end

