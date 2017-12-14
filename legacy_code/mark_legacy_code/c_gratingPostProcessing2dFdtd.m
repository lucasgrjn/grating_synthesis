classdef c_gratingPostProcessing2dFdtd
    %C_POSTPROCESSING2DFDTD Processes 2d fdtd outputs from grating
    %   Detailed explanation goes here
    
    properties
        inputs
        constants
        data   
        calc
    end
    
    methods
        function obj = c_gratingPostProcessing2dFdtd(varargin)
            if nargin == 0
                error('no inputs')
            end            
            
            domainFields = {'P'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case domainFields
                        eval(['obj.inputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            
            obj.constants.Zo = 299792458 * 4e-7*pi; % Impedance of free space = c * mu_not ~ 377ohms
        end
        
        function obj = loadData(obj)
            % load all relevant data from fdtd output files
            
            % dielectric profile
            diel = load('diel');
            obj.data.x1 = ((1:size(diel,2))-0.5)*obj.inputs.dx; 
            obj.data.y1 = ((1:size(diel,1))-0.5)*obj.inputs.dy; 
            obj.data.refIndex = sqrt(diel); % refractive index
            
            % flux and overlaps
            obj.flux_ez = load('flux_ez');
            obj.ovlpa_ez = load('ovlpa_ez');
            obj.ovlpb_ez = load('ovlpb_ez');
            
            % pull out lambda vector
            obj.data.lambda = flux_ez(:,1);
            
            % up and down flux vectors
            obj.data.up_flux = flux_ez(:,1); 
            obj.data.down_flux = -flux(:,4); % why minus? 
            
            % load... not sure
            fz03 = load('fz03'); fx03 = load('fx03'); DIRb = 1; % Downward (DIR flips the fields properly..)
            fz04 = load('fz04'); fx04 = load('fx04'); DIRt = -1; % upward
            obj.data.fz_top =      complex(fz04(:,1:2:end-1),fz04(:,2:2:end)); 
            obj.data.fx_top = DIRt*complex(fx04(:,1:2:end-1),fx04(:,2:2:end));       % [MP0] This is weird.. should be negative for bottom if anything..
            obj.data.fz_bot =      complex(fz03(:,1:2:end-1),fz03(:,2:2:end));
            obj.data.fx_bot = DIRb*complex(fx03(:,1:2:end-1),fx03(:,2:2:end));
            
            obj.data.idxLam = round( length(obj.inputs.wavelengthVec)/2); %index of center wavelength
            
            
        end
        
        function obj = plotIndex(obj)
        end
        
        function obj = calcRTFlux(obj)
            % calculate the guided reflection and transmission as and the
            % flux through all planes
            
            % Reflection and transmission in the guided mode
            f1 = complex(obj.data.ovlpa_ez(:,2), obj.data.ovlpa_ez(:,3));
            f2 = complex(obj.data.ovlpa_ez(:,4), obj.data.ovlpa_ez(:,5));
            b1 = complex(obj.data.ovlpb_ez(:,2), obj.data.ovlpb_ez(:,3));
            obj.calc.T  = abs(f2./f1).^2;
            obj.calc.R  = abs(b1./f1).^2;
            leftrad_flux  = -( flux_ez(:,2) - (conj(f1).*f1) + (conj(b1).*b1) );
            rightrad_flux =  ( flux_ez(:,3) - (conj(f2).*f2) );
            
            % up/down flux curves and directionality
            obj.calc.up_flux_norm   = obj.data.up_flux      ./ (conj(f1).*f1);
            obj.calc.down_flux_norm = obj.data.down_flux    ./ (conj(f1).*f1);
            obj.calc.directionality = obj.calc.up_flux_norm ./ obj.calc.down_flux_norm;
            obj.calc.left_flux_norm = leftrad_flux ./ (conj(f1).*f1);
            obj.calc.right_flux_norm= rightrad_flux./ (conj(f1).*f1);
        end
        
        function plotEzSnapshot(obj, whichSnapshot)
            % load and plot whichSnapshot
            % example: whichSnapshot = 'ez003'
        end
        
        
    end
    
end

