classdef c_fdtdGrating < fdtdSim
    %C_FDTDGRATING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % parameters
        inputs
        process
    end
    
    methods
        function obj = c_fdtdGrating(fdtdInputs,varargin)
            if nargin == 0
                error('no inputs')
            end
            obj = obj@fdtdSim(fdtdInputs{:});
            
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
            
            obj.P = obj.inputs.P;
        end
        
        function obj = createGrating(obj)
            % calculate dimensions
            xStartWg = 1.0;
            maxPeriods = 25;
            xmax = xStartWg; %start with 1 micron extra waveguide
            for II=1:length(obj.P.dims)
                if II<maxPeriods
                    xmax = xmax+obj.P.dims{II}.period;
                end
            end
            ymax = obj.P.thickness.sum;
            % need to set the domain in the fdtd properties
            obj.domain.domain = [xmax ymax];
            
            % layer stack
            layerStack = [obj.P.thickness.air obj.P.thickness.sio_bot obj.P.thickness.csi  obj.P.thickness.sin obj.P.thickness.psi  obj.P.thickness.sio_top];
            layers = [-1, sum(layerStack(1:1)), sum(layerStack(1:2)), sum(layerStack(1:3)),...
                sum(layerStack(1:4)),sum(layerStack(1:5))];
            indexStack = [obj.P.index.air, obj.P.index.sio_bot, obj.P.index.csi, obj.P.index.psi,...
                obj.P.index.sin,obj.P.index.sio_top];
            
            % layers
            obj = obj.addLayers('numLayers',length(indexStack([1:2 5:6])), ...
                'minY', layers([1:2 4 5]), ...
                'indices', indexStack([1:2 5:6]));
            
            xPos = xStartWg;
            
            % add starting c-Si waveguide
            obj = obj.addStructure('Rectangle', ...
                'anchorPoint',[-100 layers(3)], ...
                'structureDimensions',[100+xStartWg obj.P.thickness.csi], ...
                'index', obj.P.index.csi);
            
            % start adding unit cells
            for II=1:length(obj.P.dims)
                if II<maxPeriods %max number of periods
                    % c-Si
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint', [xPos+obj.P.dims{II}.period-obj.P.dims{II}.cSiWidth layers(3)], ...
                        'structureDimensions', [obj.P.dims{II}.cSiWidth obj.P.thickness.csi], ...
                        'index', obj.P.index.csi);
                    % conformal sin
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint', [xPos+obj.P.dims{II}.period-obj.P.dims{II}.offset-obj.P.dims{II}.pSiWidth-layerStack(4) layers(4)], ...
                        'structureDimensions', [obj.P.dims{II}.pSiWidth+2*layerStack(4) obj.P.thickness.psi+layerStack(4)], ...
                        'index', obj.P.index.sin);
                    % p-Si
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint', [xPos+obj.P.dims{II}.period-obj.P.dims{II}.offset-obj.P.dims{II}.pSiWidth layers(4)], ...
                        'structureDimensions', [obj.P.dims{II}.pSiWidth obj.P.thickness.psi], ...
                        'index', obj.P.index.psi);
                    
                    
                    
                    xPos = xPos+obj.P.dims{II}.period;
                end
            end
            
            % set values for source insertion
            sourceX = 0.1;
            sourceY1 = layers(3)-0.5;
            sourceY2 = layers(3)+obj.P.thickness.csi+0.5;
            
            
            %source from left
            obj = obj.addSource('timeType',obj.P.timeType, ...
                'spatialType','modeSolver',...
                'wavelength',obj.P.l0, ...
                'pulseWidth', obj.P.pulseWidth, ...
                'pulseDelay', obj.P.pulseDelay, ...
                'extents',[sourceX sourceX sourceY1 sourceY2], ...
                'modeNumber', 1, ...
                'nEffGuess', 2.3, ...
                'propDirection', 0);
            
            %observation planes
            %left of simulation space, after source
            obj = obj.addObsPlane('extents',[sourceX+0.1 sourceX+0.1 sourceY1 sourceY2], ...
                'typeOfUse','modeOverlap',...
                'sourceToOverlap',1);
            %right of simulation space
            obj = obj.addObsPlane('extents',[xmax-.1 xmax-.1 sourceY1 sourceY2], ...
                'typeOfUse','modeOverlap',...
                'sourceToOverlap',1);
            %bottom of simulation space
            obj = obj.addObsPlane('extents',[0 xmax .1 .1], ...
                'typeOfUse','totalPower');
            %top of simulation space
            obj = obj.addObsPlane('extents',[0 xmax ymax-0.1 ymax-0.1], ...
                'typeOfUse','totalPower');
            
        end
        
        function obj = getParams(obj)
            % get thickness and refractive index information for all layers
            % extract some inputs
            obj.P.d = obj.domain.discretization(1);
            obj.P.lambda0 = obj.inputs.P.l0;
            
            % set layer thicknesses
            obj.P.thickness.csi = roundToGrid(f_dims_12soi('bodySi','t'),obj.P.d); %crystalline Si
            obj.P.thickness.psi = roundToGrid(f_dims_12soi('polySi','t'),obj.P.d); %p-Si
            obj.P.thickness.sio_bot = roundToGrid(f_dims_12soi('STI','t'),obj.P.d); %buried oxide
            obj.P.thickness.sio_top = 1.0; %backend oxide
            obj.P.thickness.sin = roundToGrid(f_dims_12soi('SiN','t'),obj.P.d); %silicon nitride liner
            obj.P.thickness.air = 2.0;
            obj.P.thickness.pml = 0.5; %pml thickness
            obj.P.thickness.vec = [obj.P.thickness.air, obj.P.thickness.sio_bot, obj.P.thickness.csi, obj.P.thickness.psi, obj.P.thickness.sin, obj.P.thickness.sio_top]; %vector of layer thicknesses
            obj.P.thickness.sum = sum(obj.P.thickness.vec); % total thickness
            
            % set refractive indexes
            obj.P.index.sio_bot = index_IBM12SOI45nm_fits(obj.P.lambda0,'STI'); %buried oxide
            obj.P.index.sio_top = index_IBM12SOI45nm_fits(obj.P.lambda0,'PSG'); %backend oxide
            obj.P.index.sin = index_IBM12SOI45nm_fits(obj.P.lambda0,'SiN'); %silicon nitride liner
            obj.P.index.csi = index_IBM12SOI45nm_fits(obj.P.lambda0,'bodySi'); %crystalline silicon
            obj.P.index.psi = index_IBM12SOI45nm_fits(obj.P.lambda0,'polySi'); %polysilicon
            obj.P.index.air = 1.0;
        end
        
        function obj = processData(obj,varargin) %this is terrible! fix!
            processingFields = {'OPTS'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case processingFields
                        eval(['obj.process.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            if ~isfield(obj.process.OPTS, 'fibThetaScan')
                obj.process.OPTS.fibThetaScan = 5:2:30;
            end
            if ~isfield(obj.process.OPTS, 'xOffsetScan')
                obj.process.OPTS.xOffsetScan = 0:0.1:10;
            end
            [fig_file,fig_ext,fignum,directionality,Qgp] = getParams3mod_1280nm_down_markEdit_v4('dx',obj.domain.discretization(1),'dy',obj.domain.discretization(2),'wavelengthVec',obj.domain.wavelengthSpectrum(1): obj.domain.wavelengthSpectrum(3): obj.domain.wavelengthSpectrum(2),'xmax',obj.domain.domain(1),'MFD',10,'direction','down','OPTS',obj.process.OPTS, ...
                'fibThetaScan',obj.process.OPTS.fibThetaScan,'xOffsetScan',obj.process.OPTS.xOffsetScan);
            if ~strcmp(obj.process.OPTS.switches.justLoadData,'yes')
                name = obj.domain.nameDir;
                nameBase = obj.domain.nameDirBase;
                %                 FID = fopen(strcat(name,'_data_.txt'),'w');
                
                mkdir(name)
                if Qgp.OPTS.FSAVEFIGS
                    for ii=find(Qgp.OPTS.OUTPLOTS.*(1:length(Qgp.OPTS.OUTPLOTS)))
                        movefile(strcat(fig_file,num2str(ii),'.',fig_ext),name);
                    end
                end
                obj.moveFiles;
                save([nameBase '_gratSimObj'],'obj')
            end
            obj.postProcessingResults = Qgp;
        end
    end
    
end

