classdef c_twoLevelGratingCell < c_gratingCell
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
        
    end     % end properties
    
    
    methods
        
        function obj = c_twoLevelGratingCell(varargin)
            % Constructor
            
            % call grating cell constructor
            obj = obj@c_gratingCell(varargin{:});
            
        end     % end function addLayer()
        
        
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
            
        end     % end function twoLevelBuilder
   
        
    end     % end methods
    
end     % end class





























