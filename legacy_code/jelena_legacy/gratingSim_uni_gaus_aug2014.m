classdef gratingSim_uni_gaus_aug2014 < fdtdSim
    %this is the help file
    
    properties
        gratingDesign
        gratingDims
        gratingInputs
        gratingDielectrics
        processingInputs
    end
    
    methods
        function obj = gratingSim_uni_gaus_aug2014(fdtdInputs,varargin)
            if nargin == 0
                error('no inputs')
            end
            obj = obj@fdtdSim(fdtdInputs{:});
            
            domainFields = {'P'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case domainFields
                        eval(['obj.gratingInputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end            
        end
        
        function obj = createGrating(obj,varargin)
            % builds the grating simulation (structures, sources, observation planes, etc.)
            %% parse inputs
            gratingFields = {'process','whichDesign','OPTS','varLoops'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case gratingFields
                        eval(['obj.gratingDielectrics.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            %% get design
            if ~isfield(obj.gratingDesign,'whichDesign')
                error('No grating design specified')
            end
            obj = obj.getGratingDesign;
            %% build dielectric profile
            switch obj.gratingDesign.whichDesign
                %case 'unidirDown1280EOS18_TEM' %build dielectric cross section for IBM
                %case 'unidirDown1200EOS20' %build dielectric cross section for IBM
                %case 'uniform1280_EOS18_20_designDims' %build dielectric cross section for IBM
                case 'unidirDown1280EOS24Gaussian' %build dielectric cross section for IBM
                %case 'uniform1280_EOS18_20_TEMDims' %build dielectric cross section for IBM
                %case 'uniform1200EOS22backward' %build dielectric cross section for IBM    
                %case 'uniform1200EOS22backwardART' %build dielectric cross section for IBM   
                %case 'uniform1280EOS22backwardART' %build dielectric cross section for IBM  
                %case 'uniform1200EOS22forward' %build dielectric cross section for IBM
                %case 'uniform1200EOS22forwardART' %build dielectric cross section for IBM
                %case 'uniform1280EOS22forward' %build dielectric cross section for IBM
                %case 'uniform1280EOS22forwardART' %build dielectric cross section for IBM
                    %grab a few variable values
                    if ~isfield(obj.gratingDesign.OPTS,'psiAbsorption')
                        psiAbsorption = 0; 
                    else
                        psiAbsorption = obj.gratingDesign.OPTS.psiAbsorption;
                    end
                    if ~isfield(obj.gratingDesign.OPTS,'addImperfections')
                        addImperfections = 'yes';
                    else
                        addImperfections = obj.gratingDesign.OPTS.addImperfections;
                    end                    
                    
                    obj = obj.varLoops;
                    
                    TL = obj.gratingDims.TL;

                    xmax = obj.domain.domain(1);
                    ymax = obj.domain.domain(2);
                    
                    dx = obj.domain.discretization(1); % [MP] temporarily define discretization to use below.
                    
                    
                    layerstack = [TL.air_thickness, TL.sio_bottom, TL.si_thickness, ...
                        TL.psi_thickness, TL.sin_thickness];
                    layers = [-1, sum(layerstack(1:1)), sum(layerstack(1:2)), sum(layerstack(1:3)),...
                        sum(layerstack(1:4)),sum(layerstack(1:5))];
                    indexstack = [TL.air_index, TL.sio_index_1550, TL.si_index_1550, TL.psi_index_1550,...
                        TL.sin_index_1550,TL.sio_index_1550];
                    
                    
                    obj = obj.addLayers('numLayers',length(indexstack), ...
                        'minY', layers, ...
                        'indices', indexstack);
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 layers(4)], ...
                        'structureDimensions',[2*xmax TL.sin_thickness+TL.psi_thickness], ...
                        'index',TL.sio_index_1550);
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 layers(4)], ...
                        'structureDimensions',[2*xmax TL.sin_thickness], ...
                        'index',TL.sin_index_1550);
                    
                    if strcmp(addImperfections,'no')
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[-5 layers(4)], ...
                            'structureDimensions',[2*xmax TL.sin2ImpThickness], ...
                            'index',TL.sin2_index_1550);
                    end
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 -5], ...
                        'structureDimensions',[2*xmax 5+TL.air_thickness], ...
                        'index',TL.air_index);
                    
                   
                    %ypos = 5.0; %start grating ARtooth and unit cell teeth
                    ypos = 2.0; %start grating ARtooth and unit cell teeth
%                     obj = obj.addStructure('Rectangle', ...
%                         'anchorPoint',[ypos+TL.aroff+TL.argap layers(3)], ...
%                         'structureDimensions',[TL.argap TL.si_thickness], ...
%                         'index',TL.sio_index_1550);
                    % adding non uniform periods in C-Si layer for gaussian Shape

                    
                    
%%%%%%%%%%%%%%% CHANGES STARTING HERE!!!!!!!!
                    if TL.barpoly(1) < .5*TL.GPbar(1)
                        inverteddesign = 1;
                    else
                        inverteddesign = 0;
                    end
                    
                    % si gap objects (first is different then loop for rest)
                    y1 = round( ypos/dx ) * dx;
                    y2 = round( (ypos+TL.GPgap(1)) / dx ) * dx;
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[y1 layers(3)], ...
                        'structureDimensions',[(y2-y1) TL.si_thickness], ...
                        'index',TL.sio_index_1550);
                    for j = 2:length(TL.GPgap)
                        y1 = round( (ypos+sum([TL.GPgap(1:j-1) TL.GPbar(1:j-1)])) / dx ) * dx;
                        y2 = round( (ypos+sum([TL.GPgap(1:j)   TL.GPbar(1:j-1)])) / dx ) * dx;
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[y1 layers(3)], ...
                            'structureDimensions',[(y2-y1) TL.si_thickness], ...
                            'index',TL.sio_index_1550);
                    end
                    
%                     % si bar objects (first then loop for rest)
%                     obj = obj.addStructure('Rectangle', ...
%                         'anchorPoint',[ypos+TL.GPgap(1) layers(3)], ...
%                         'structureDimensions',[TL.GPbar(1) TL.si_thickness], ...
%                         'index',TL.si_index_1550);
%                     for j = 2:length(TL.GPbar)
%                         obj = obj.addStructure('Rectangle', ...
%                         'anchorPoint',[ypos+sum([TL.GPgap(1:j) TL.GPbar(1:j-1)]) layers(3)], ...
%                         'structureDimensions',[TL.GPbar(j) TL.si_thickness], ...
%                         'index',TL.si_index_1550);
%                     end
                    
                    % input p-si waveguide SiN ribbon ontop of c-si
                    if inverteddesign == 0
                        y1 = -round(1.0/dx)*dx; % Start at -1um.
                        y2 = round( (ypos-TL.polyoffset+TL.sin_thickness)/dx ) * dx;
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[-1 layers(4)], ...
                            'structureDimensions',[(y2-y1) TL.psi_thickness+TL.sin_thickness], ...
                            'index',TL.sin_index_1550);
%                     else
%                         obj = obj.addStructure('Rectangle', ...
%                             'anchorPoint',[0 layers(4)], ...
%                             'structureDimensions',[TL.sin_thickness+ypos-TL.offset(1) TL.psi_thickness+TL.sin_thickness], ...
%                             'index',TL.sin_index_1550);
                    end
                
                    % p-si bar SiN objects with sin ribbon
                    for j = 1:length(TL.barpoly)
                        y1 = round((ypos + sum([TL.GPgap(1:j) TL.GPbar(1:j)]) - TL.offset(j) - TL.barpoly(j) - TL.sin_thickness)/dx) * dx;
                        y2 = round((ypos + sum([TL.GPgap(1:j) TL.GPbar(1:j)]) - TL.offset(j)                 + TL.sin_thickness)/dx) * dx;
                        obj = obj.addStructure('Rectangle', ...
                       'anchorPoint',[y1 layers(4)], ...
                       'structureDimensions',[(y2-y1) TL.psi_thickness+TL.sin_thickness], ...
                       'index',TL.sin_index_1550);
                    end
                    % p-si bar poly object
                    for j = 1:length(TL.barpoly)
                        y1 = round((ypos + sum([TL.GPgap(1:j) TL.GPbar(1:j)]) - TL.offset(j) - TL.barpoly(j))/dx) * dx;
                        y2 = round((ypos + sum([TL.GPgap(1:j) TL.GPbar(1:j)]) - TL.offset(j))/dx) * dx;
                        obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[y1 layers(4)], ...
                        'structureDimensions',[(y2-y1) TL.psi_thickness], ...
                        'index',TL.psi_index_1550,...
                        'absorption',psiAbsorption);
                    end
                    
%                     % % %
%                     % tooth before start
%                     y1 = round((ypos - TL.offset(1) - TL.sin_thickness)/dx) * dx;
%                     y2 = round((ypos - TL.offset(1)                 + TL.sin_thickness)/dx) * dx;
%                     obj = obj.addStructure('Rectangle', ...
%                         'anchorPoint',[y1 layers(4)], ...
%                         'structureDimensions',[(y2-y1) TL.psi_thickness+TL.sin_thickness], ...
%                         'index',TL.sin_index_1550);
% %                     y1 = round((ypos - TL.offset(1) - TL.barpoly(1))/dx) * dx;
% %                     y2 = round((ypos - TL.offset(1))/dx) * dx;
% %                     obj = obj.addStructure('Rectangle', ...
% %                         'anchorPoint',[y1 layers(4)], ...
% %                         'structureDimensions',[(y2-y1) TL.psi_thickness], ...
% %                         'index',TL.psi_index_1550,...
% %                         'absorption',psiAbsorption);
%                     % % %
                    
                    % input p-si waveguide ontop of c-si
                    if inverteddesign == 0
                        y1 = -round(1.0/dx)*dx; % Start at -1um.
                        y2 = round( (ypos-TL.polyoffset)/dx ) * dx;
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[y1 layers(4)], ...
                            'structureDimensions',[(y2-y1) TL.psi_thickness], ...
                            'index',TL.psi_index_1550,...
                            'absorption',psiAbsorption);
                    end
                    
%%%%%%%%%%%%%%% CHANGES ENDING HERE!!!!!!!!

                   %----------------------------------------------------------
 
                   %obj = obj.addStructure('Rectangle', ...
                        %'anchorPoint',[ypos-TL.GPbar-TL.GPgap layers(3)], ...
                        %'structureDimensions',[TL.GPgap TL.si_thickness], ...
                        %'index',TL.sio_index_1550);
                    %obj = obj.addStructure('Rectangle', ...
                        %'anchorPoint',[ypos-TL.aroff-TL.argap layers(3)], ...
                        %'structureDimensions',[.1 TL.si_thickness], ...
                        %'index',TL.sio_index_1550);
                    
                    for ii = 1:TL.periods
                        period = TL.gap+TL.bar;
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos layers(3)], ...
                            'structureDimensions',[TL.gap TL.si_thickness], ...
                            'index',TL.sio_index_1550);
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos-TL.barpoly-TL.offset+period layers(4)], ...
                            'structureDimensions',[TL.barpoly TL.psi_thickness], ...
                            'index',TL.psi_index_1550,...
                            'absorption',psiAbsorption);
                        
                        if strcmp(addImperfections,'yes') %if simulating imperfections, change where main sin height is
                            sinOffset = TL.sin2ImpThickness;
                        else
                            sinOffset = 0.0;
                        end
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos-TL.barpoly-TL.offset-TL.sin_thickness+period layers(4)+sinOffset], ...
                            'structureDimensions',[TL.sin_thickness TL.psi_thickness+TL.sin_thickness-sinOffset], ...
                            'index',TL.sin_index_1550);
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos-TL.barpoly-TL.offset+period layers(4)+TL.psi_thickness], ...
                            'structureDimensions',[TL.barpoly TL.sin_thickness], ...
                            'index',TL.sin_index_1550);
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos-TL.offset+period layers(4)+sinOffset], ...
                            'structureDimensions',[TL.sin_thickness TL.psi_thickness+TL.sin_thickness-sinOffset], ...
                            'index',TL.sin_index_1550);
                        
                        if strcmp(addImperfections,'yes')
                            
                            %add dip in pSi from TEM
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.barpoly-TL.offset+period+TL.barpoly-(TL.bar-TL.offset)-TL.pSiImpW layers(4)-TL.pSiImpH], ...
                                'structureDimensions',[TL.pSiImpW TL.pSiImpH], ...
                                'index',TL.psi_index_1550,...
                                'absorption',psiAbsorption);
                            %add dip in SiN by bodySi
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos layers(4)-TL.sinImpH], ...
                                'structureDimensions',[TL.sinImpW TL.sinImpH], ...
                                'index',TL.sin2_index_1550);
                            %add 3 SiN liner sections
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.barpoly-TL.offset-TL.sin2ImpThickness+period layers(4)], ...
                                'structureDimensions',[TL.sin2ImpThickness TL.psi_thickness+TL.sin2ImpThickness], ...
                                'index',TL.sin2_index_1550);
                            
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.barpoly-TL.offset+period layers(4)+TL.psi_thickness], ...
                                'structureDimensions',[TL.barpoly TL.sin2ImpThickness], ...
                                'index',TL.sin2_index_1550);
                            
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.offset+period layers(4)], ...
                                'structureDimensions',[TL.sin2ImpThickness TL.psi_thickness+TL.sin2ImpThickness], ...
                                'index',TL.sin2_index_1550);
                            
                            %add oxide spacers
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.barpoly-TL.offset+period-TL.widthSpacer layers(4)], ...
                                'structureDimensions',[TL.widthSpacer TL.psi_thickness], ...
                                'index',TL.sio_index_1550,...
                                'absorption',psiAbsorption);
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos-TL.barpoly-TL.offset+period+TL.barpoly layers(4)], ...
                                'structureDimensions',[TL.widthSpacer TL.psi_thickness], ...
                                'index',TL.sio_index_1550,...
                                'absorption',psiAbsorption);
                            
                        end
                        
                        ypos = ypos+TL.gap+TL.bar;
                    end
                    
                    %source from left
                    obj = obj.addSource('timeType',obj.gratingInputs.P.timeType, ...
                        'spatialType','modeSolver',...
                        'wavelength',obj.gratingInputs.P.l0, ...
                        'pulseWidth', 10, ...
                        'pulseDelay', 3, ...
                        'extents',[.1 .1 0 ymax], ...
                        'modeNumber', 1, ...
                        'nEffGuess', 1.7, ...
                        'propDirection', 0);
                    
                    %observation planes
                    %left of simulation space, after source
                    obj = obj.addObsPlane('extents',[.2 .2 0 ymax], ...
                        'typeOfUse','modeOverlap',...
                        'sourceToOverlap',1);
                    %right of simulation space
                    obj = obj.addObsPlane('extents',[xmax-.1 xmax-.1 0 ymax], ...
                        'typeOfUse','modeOverlap',...
                        'sourceToOverlap',1);
                    %bottom of simulation space
                    obj = obj.addObsPlane('extents',[0 xmax .1 .1], ...
                        'typeOfUse','totalPower');
                    %top of simulation space
                    obj = obj.addObsPlane('extents',[0 xmax ymax ymax], ...
                        'typeOfUse','totalPower');
                
                case {'uniform1280_EOS18_20_designDims','uniform1280_EOS18_TEMDims'}    
                    TL = obj.gratingDims.TL;
                    if ~isfield(obj.gratingDesign.OPTS,'psiAbsorption')
                        psiAbsorption = 0;
                    else
                        psiAbsorption = obj.gratingDesign.OPTS.psiAbsorption;
                    end
                    if ~isfield(obj.gratingDesign.OPTS,'addImperfections')
                        addImperfections = 'yes';
                    else
                        addImperfections = obj.gratingDesign.OPTS.addImperfections;
                    end
                    if isfield(obj.gratingDielectrics,'varLoops')
                        if isfield(obj.gratingDesign.varLoops,'TL')
                            varFNTL = fieldnames(obj.gratingDesign.varLoops.TL);
                            switch obj.gratingDesign.varLoops.varType
                                case 'percent'
                                    for ii = 1:length(varFNTL)
                                        eval(['TL.' varFNTL(ii) '=TL.' varFNTL(ii) '*(1+obj.gratingDesign.varLoops.TL.' varFNTL(ii) ')'])
                                    end
                                case 'absolute'
                                    for ii = 1:length(varFNTL)
                                        eval(['TL.' varFNTL(ii) '= obj.gratingDesign.varLoops.TL.' varFNTL(ii) ')'])
                                    end
                                otherwise
                                    error('Specify type of variation. percent or absolute')
                            end
                        end
                        
                    end
                    
                    xmax = obj.domain.domain(1);
                    ymax = obj.domain.domain(2);
                    
                    layerstack = [TL.air_thickness, TL.sio_bottom, TL.si_thickness, ...
                        TL.psi_thickness, TL.sin_thickness];
                    layers = [-1, sum(layerstack(1:1)), sum(layerstack(1:2)), sum(layerstack(1:3)),...
                        sum(layerstack(1:4)),sum(layerstack(1:5))];
                    indexstack = [TL.air, TL.sio_index_1550, TL.si_index_1550, TL.psi_index_1550,...
                        TL.sin_index_1550,TL.sio_index_1550];
                    
                    %layers
                    obj = obj.addLayers('numLayers',length(indexstack), ...
                        'minY', layers, ...
                        'indices', indexstack);
                    %remove unneeded (for now) layers
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 layers(4)], ...
                        'structureDimensions',[50 TL.sin_thickness+TL.psi_thickness], ...
                        'index',TL.sio_index_1550);
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 layers(4)], ...
                        'structureDimensions',[50 TL.sin_thickness], ...
                        'index',TL.sin_index_1550);
                    
                    obj = obj.addStructure('Rectangle', ...
                        'anchorPoint',[-5 -5], ...
                        'structureDimensions',[50 5+TL.air_thickness], ...
                        'index',TL.air);
                    
                    ypos = 2.0;
                    for ii = 1:TL.periods
                        period = TL.gap+TL.bar;
                        
                        obj = obj.addStructure('Rectangle', ...
                            'anchorPoint',[ypos layers(3)], ...
                            'structureDimensions',[TL.gap TL.si_thickness], ...
                            'index',TL.sio_index_1550);
                        if strcmp(addImperfections,'yes')
                            
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos layers(4)-TL.sinImpH], ...
                                'structureDimensions',[TL.sinImpW TL.sinImpH], ...
                                'index',TL.sin_index_1550);
                            
                            obj = obj.addStructure('Rectangle', ...
                                'anchorPoint',[ypos+TL.gap-TL.sinImpW layers(4)-TL.sinImpH], ...
                                'structureDimensions',[TL.sinImpW TL.sinImpH], ...
                                'index',TL.sin_index_1550);
                        end
                        ypos = ypos+TL.gap+TL.bar;
                    end
                    
                    %source from left
                    obj = obj.addSource('timeType',obj.gratingInputs.P.timeType, ...
                        'spatialType','modeSolver',...
                        'wavelength',obj.gratingInputs.P.l0, ...
                        'pulseWidth', 10, ...
                        'pulseDelay', 3, ...
                        'extents',[.1 .1 0 ymax], ...
                        'modeNumber', 1, ...
                        'nEffGuess', 1.7, ...
                        'propDirection', 0);
                    
                    %observation planes
                    %left of simulation space, after source
                    obj = obj.addObsPlane('extents',[.2 .2 0 ymax], ...
                        'typeOfUse','modeOverlap',...
                        'sourceToOverlap',1);
                    %right of simulation space
                    obj = obj.addObsPlane('extents',[xmax-.1 xmax-.1 0 ymax], ...
                        'typeOfUse','modeOverlap',...
                        'sourceToOverlap',1);
                    %bottom of simulation space
                    obj = obj.addObsPlane('extents',[0 xmax 0 0], ...
                        'typeOfUse','totalPower');
                    %top of simulation space
                    obj = obj.addObsPlane('extents',[0 xmax ymax ymax], ...
                        'typeOfUse','totalPower');
                    
                otherwise
                    error('No valid process specified')
            end
            
        end
        
        function obj = runGratingSimulation(obj,varargin)
            
            gratingFields = {'process','whichDesign','OPTS','varLoops'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case gratingFields
                        eval(['obj.gratingDesign.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            
            if isfield(obj.gratingDesign,'varLoops')
                obj.gratingDesign.varLoopsBase = obj.gratingDesign.varLoops;
                obj.domain.nameddBase = obj.domain.namedd;
                obj.domain.nameDirBase = ['LOOPS_' obj.domain.nameDir]; %create unchanging root directory for looped data
                fieldNTL = fieldnames(obj.gratingDesign.varLoops.TL);               
                for ii = 1:length(fieldNTL)
                    T(ii) = eval(['length(obj.gratingDesign.varLoops.TL.' fieldNTL{ii} ')']);
                end               
                obj = obj.generateLoops(T,length(fieldNTL));
            else
                obj.domain.nameddBase = obj.domain.namedd;
                obj.domain.nameDirBase = obj.domain.nameDir; %create unchanging root directory for looped data                
                obj = obj.createGrating('process',obj.gratingDesign.process,'whichDesign',obj.gratingDesign.whichDesign,...
                    'OPTS',obj.gratingDesign.OPTS);
                obj.runSimulation
                obj.processData('OPTS',obj.gratingDesign.OPTS)                
            end
            
        end
        
        function generateReport %intended to output a smart text file and/or data file that contains all info you'd want from variants etc
        end
        
        function obj = processData(obj,varargin)
            processingFields = {'OPTS'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case processingFields
                        eval(['obj.processingInputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end            
            %[fig_file,fig_ext,fignum,directionality,Qgp] = getParams3mod_1280nm_down_markEdit_v4('dx',obj.domain.discretization(1),'dy',obj.domain.discretization(2),'wavelengthVec',obj.domain.wavelengthSpectrum(1): obj.domain.wavelengthSpectrum(3): obj.domain.wavelengthSpectrum(2),'xmax',obj.domain.domain(1),'MFD',5,'direction','down','OPTS',obj.processingInputs.OPTS, ...
                %'fibThetaScan',5:0.5:30,'xOffsetScan',0:0.1:6);
                mfd = 5;
                offsetscan = 0:0.1:10;
            [fig_file,fig_ext,fignum,directionality,Qgp] = getParams3mod_1280nm_down_markEdit_v4RK('dx',obj.domain.discretization(1),'dy',obj.domain.discretization(2),'wavelengthVec',obj.domain.wavelengthSpectrum(1): obj.domain.wavelengthSpectrum(3): obj.domain.wavelengthSpectrum(2),'xmax',obj.domain.domain(1),'MFD',mfd,'direction','down','OPTS',obj.processingInputs.OPTS, ...
                'fibThetaScan',-40:1:40,'xOffsetScan',offsetscan);
                if ~strcmp(obj.processingInputs.OPTS.switches.justLoadData,'yes')
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
                    save(['gratSimObj_' nameBase],'obj')
                end
            
            
        end
        
        function obj = generateLoops(obj,T,N,INDEXVEC,NMAX,fieldNTL)
            % recursive function to generate arbitrary number of simulation
            % for loops
            if nargin == 3 || isempty(INDEXVEC)
                INDEXVEC = ones(1,N); %initialize index vector
                NMAX = N;
                fieldNTL = fieldnames(obj.gratingDesign.varLoops.TL);
            end
            
            if N < 1
                error(sprintf('ITERATE called with N < 0 (N = %d)',N));
            end
            
            for ii=1:T(length(T)-N+1)
                INDEXVEC(length(INDEXVEC)-N+1) = ii;
                switch N
                    case 1
                        fprintf(['INDEXVEC: ', repmat('%d',1,length(INDEXVEC)),'\n'],INDEXVEC)
                        for jj=1:NMAX
                            fprintf([fieldNTL{jj} '=' num2str(eval(['obj.gratingDesign.varLoops.TL.' fieldNTL{jj} '(INDEXVEC(jj))'])) '\n'])
                            eval(['obj.gratingDesign.varLoops.TL.' fieldNTL{jj} '_ITER = obj.gratingDesign.varLoops.TL.' fieldNTL{jj} '(INDEXVEC(jj))']);
                        end
                        obj.domain.namedd = strcat(obj.domain.nameddBase,'\',strrep(num2str(INDEXVEC),' ','')); % change property to make new file for each iteration
                        obj.domain.nameDir = strcat(obj.domain.nameDirBase,'\',strrep(num2str(INDEXVEC),' ','')); % change property to make new folder for each iteration
                        obj = obj.createGrating('process',obj.gratingDesign.process,'whichDesign',obj.gratingDesign.whichDesign,...
                            'OPTS',obj.gratingDesign.OPTS);
                        obj.runSimulation
                        obj.processData('OPTS',obj.gratingDesign.OPTS)
                        obj.gratingDims = [];
                        obj.sources = [];
                        obj.structures = [];
                        obj.obsPlanes = [];
                        
                    otherwise
                        obj = obj.generateLoops(T,N-1,INDEXVEC,NMAX,fieldNTL);
                end
            end
            %write loop info to a text file
            nameBase = obj.domain.nameDirBase;
            fid = fopen([nameBase '\loopInfo.txt'],'wt');
            fieldNTL = fieldnames(obj.gratingDesign.varLoopsBase.TL);
            tempStr = ['variation type = ' obj.gratingDesign.varLoops.varType '\n'];
            for ii=1:length(fieldNTL)
                tempStr = [tempStr 'index ' num2str(ii) ' = ' fieldNTL{ii} ' and has value [' num2str(eval(['obj.gratingDesign.varLoops.TL.' fieldNTL{ii}]),'%g\x20') '] \n'];
            end
            fprintf(fid,tempStr);
            fclose(fid);
        end
        
        function obj = getGratingDesign(obj)
            % modifies obj.gratingDims.TL property to return layer params            
            %Input:
            %   uses obj.gratingDesign.whichDesign as switch. Valid values for obj.gratingDesign.whichDesign:
            %           unidirDown1280EOS18_20
            %           unidirDown1200EOS20
            %           unidirUp1200EOS20
            %           unidirUp1280EOS18_20
            %           uniform1200Stevan_EOS16_20
            %           uniform1280_EOS18_20_designDims
            %           uniform1280HighAngle_EOS20
            %           uniform1280HighAngle2_EOS20
            %           uniform1280D1L100umRem
            %           uniform1550D1L100umRem
            %           unidirDown1280EOS18_20_gdsTest
            %           unidirDown1280EOS18_TEM
            %           uniform1280_EOS18_TEMDims
            %
            % Output:
            %   obj.gratingDims.TL containing index of refraction and
            %   dimensions of layers
            %
            % reference: Kareem's grating design pdf
            
            dx = obj.domain.discretization(1);
            dy = obj.domain.discretization(2);
            
            switch obj.gratingDesign.whichDesign
                
                case 'unidirDown1280EOS24Gaussian'
                    obj.gratingDims.TL.air_index = 1;
                    obj.gratingDims.TL.si_index_1550 = 3.519;
                    obj.gratingDims.TL.psi_index_1550 = 3.59;
                    obj.gratingDims.TL.sio_index_1550 = 1.449;
                    obj.gratingDims.TL.sin_index_1550 = 1.98;
                    %layer thicknesses;
                    obj.gratingDims.TL.sio_bottom = .1415;
                    obj.gratingDims.TL.si_thickness = .07346;
                    obj.gratingDims.TL.psi_thickness = .081;
                    obj.gratingDims.TL.sin_thickness = .0723;
                    obj.gratingDims.TL.air_thickness = .8;
                    %tooth perameters;
                    obj.gratingDims.TL.gap = .39; %L insted of T
                    obj.gratingDims.TL.bar = .30;
                    obj.gratingDims.TL.offset = .15;
                    obj.gratingDims.TL.barpoly = .38;
                    obj.gratingDims.TL.aroff = .35;
                    obj.gratingDims.TL.argap = .500;
                    obj.gratingDims.TL.periods = 0;
                    
%%%%%%%%%%%%%%% CHANGES STARTING HERE!!!!!!!!
                    
%                     % best A design (without r) so far...
%                     cellstocutout = [1 2 3 4 5];
%                     % cellstocutout = [];
%                     a_final = [4110,4110,4110,490,505,520,535,565,595,715,745,745,745,745,745,745,745,745,745,745,745,745,745,745,745,745,745;]
%                     a_final(cellstocutout) = [];
%                     fill_final = [0.940000000000000,0.940000000000000,0.910000000000000,0.850000000000000,0.820000000000000,0.760000000000000,0.700000000000000,0.640000000000000,0.550000000000000,0.340000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000;]
%                     fill_final(cellstocutout) = [];
%                     r_final = 1*ones(1,length(a_final));
%                     off_final = [0.820000000000000,0.820000000000000,0.800000000000000,0.800000000000000,0.760000000000000,0.780000000000000,0.760000000000000,0.760000000000000,0.760000000000000,0.760000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000,0.780000000000000;]
%                     off_final(cellstocutout) = [];
                    
%                      % best B design so far...
%                      a_final = [625,625,625,625,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650;]
%                      fill_final = [0.800000000000000,0.800000000000000,0.800000000000000,0.700000000000000,0.600000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000;]
%                      r_final = [0.100000000000000,0.100000000000000,0.300000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000;]
%                      off_final = [0.200000000000000,0.200000000000000,0.100000000000000,0,0,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000,0.900000000000000;]

%                     % best A results from parallel sweep 1
%                      final_matrix = [4110,4110,490,490,505,520,535,565,595,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685;0.950000000000000,0.950000000000000,0.900000000000000,0.850000000000000,0.800000000000000,0.1100000000000000,0.700000000000000,0.650000000000000,0.550000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0.200000000000000,0.200000000000000,0.200000000000000,0.200000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;]
%                     cellstocutout = [];
%                     polyoffdefined = 0;
%                      cellstocutout = [1, 2, 3, 4];
%                     obj.gratingDims.TL.polyoffset = .1470;
%                     polyoffdefined = 1;
                    
%                     % best DRC A result (from parallel sweep 1) ( BEST!!! )
%                     final_matrix = [505,520,535,565,595,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685;0.800000000000000,0.1100000000000000,0.700000000000000,0.650000000000000,0.550000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;]
%                     cellstocutout = [];
%                     obj.gratingDims.TL.polyoffset = .1470;
%                     polyoffdefined = 1;
%                      
%                      % test A result from sweep 2
%                      final_matrix = [466,466,478,478,502,514,526,538,562,598,658,718,718,718,718,718,718,718,718,718,718,718,718,718,718,718,718,718,718;0.960000000000000,0.960000000000000,0.930000000000000,0.900000000000000,0.870000000000000,0.810000000000000,0.780000000000000,0.690000000000000,0.660000000000000,0.540000000000000,0.420000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000,0.330000000000000;1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,0.980000000000000,1.01000000000000,0.980000000000000,1.04000000000000,0.950000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000;0.190000000000000,0.190000000000000,0.190000000000000,0.230000000000000,0.210000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.270000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000;]
%                     cellstocutout = [];

%                      
%                     % best B results from parallel sweep 1
%                     final_matrix = [610,610,640,640,580,580,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685;0.900000000000000,0.900000000000000,0.850000000000000,0.800000000000000,0.1100000000000000,0.600000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000;0.100000000000000,0.100000000000000,0.150000000000000,0.250000000000000,0.600000000000000,0.950000000000000,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0.600000000000000,0.600000000000000,0.550000000000000,0.500000000000000,0.400000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;]
%                     cellstocutout = [];
%                     polyoffdefined = 0;

%                     % testing B from sweep 2 try 1 (81%)
%                     final_matrix = [610,610,622,634,598,586,610,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658;0.900000000000000,0.900000000000000,0.870000000000000,0.780000000000000,0.780000000000000,0.630000000000000,0.480000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000;0.0800000000000000,0.0800000000000000,0.140000000000000,0.350000000000000,0.440000000000000,0.860000000000000,1.13000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000;0.590000000000000,0.590000000000000,0.550000000000000,0.450000000000000,0.430000000000000,0.290000000000000,0.210000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000;]
%                     cellstocutout = [];
%                     polyoffdefined = 0;
% 
%                     % testing B from sweep 2 try 2 (87%)
%                     final_matrix = [610,622,622,610,610,598,646,646,646,646,646,646,646,646,646;...
%                         .9,.87,.81,.78,.69,.57,.45,.45,.45,.45,.45,.45,.45,.45,.45;...
%                         .08,.14,.29,.41,.65,.98,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01;...
%                         .59,.55,.49,.45,.35,.25,.23,.23,.23,.23,.23,.23,.23,.23,.23;]
%                       cellstocutout = [1];
%                     polyoffdefined = 0;

%                     % testing B from sweep 2 try 3
%                     final_matrix = [610,610,622,634,622,586,646,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658,658;0.900000000000000,0.900000000000000,0.870000000000000,0.780000000000000,0.1100000000000000,0.540000000000000,0.450000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000,0.420000000000000;0.0800000000000000,0.0800000000000000,0.140000000000000,0.350000000000000,0.440000000000000,1.16000000000000,1.01000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000,1.07000000000000;0.590000000000000,0.590000000000000,0.550000000000000,0.450000000000000,0.410000000000000,0.210000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000,0.230000000000000;];
%                     cellstocutout = [1, 2, 3];
%                     polyoffdefined = 0;

%                     % testing B from sweep 2 try 4
%                     final_matrix = [624.365054602184,623.858702919545,619.807889458435,613.2253111084132,569.172721194562,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001,585.31109110039001;0.880967238689548,0.877635391129931,0.850980610652998,0.807666592377981,0.517795854691331,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064,0.624414976599064;0.100000000000000,0.110000000000000,0.190000000000000,0.320000000000000,1.19000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000,0.870000000000000;0.574749694572572,0.569792139727983,0.5312807811012384,0.473056793297084,0.222251865984130,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133,0.2864131720110133;]
%                     cellstocutout = [];
%                      polyoffdefined = 0;

%                         % testing B from sweep 3 try 1
%                         final_matrix = [640,640,640,630,620,620,640,640,640,640,640,640,640,640,640,640,640,640,640,640;0.910000000000000,0.870000000000000,0.850000000000000,0.790000000000000,0.690000000000000,0.550000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000,0.470000000000000;0.0700000000000000,0.120000000000000,0.220000000000000,0.350000000000000,0.580000000000000,0.880000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000,0.940000000000000;0.600000000000000,0.560000000000000,0.530000000000000,0.480000000000000,0.370000000000000,0.280000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000;]
%                         polyoffdefined = 0;
%                         cellstocutout = [1 2];
%                         
                        %final_matrix = [640,640,640,640,640,640,640,640,640,640,640,640,640,640;0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000,0.910000000000000;0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000,0.0700000000000000;0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000,0.600000000000000;]
                        

%                     % testing B from sweep 2 try 2 (87%)
%                     final_matrix = [610,622,622,610,610,598,646,646,646,646,646,646,646,646,646;...
%                         .9,.87,.81,.78,.69,.57,.45,.45,.45,.45,.45,.45,.45,.45,.45;...
%                         .08,.14,.29,.41,.65,.98,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01;...
%                         .59,.55,.49,.45,.35,.25,.23,.23,.23,.23,.23,.23,.23,.23,.23;];
%                       cellstocutout = [];
%                     polyoffdefined = 0;
%                     

%                     % testing B from sweep 3 try 2
%                     final_matrix = [634,622,598,598,646,682,682,682,682,682,682,682,682,682,682,682,682,682;0.840000000000000,0.810000000000000,0.1100000000000000,0.570000000000000,0.450000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000,0.390000000000000;0.170000000000000,0.290000000000000,0.500000000000000,0.950000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000,1.01000000000000;0.550000000000000,0.490000000000000,0.410000000000000,0.250000000000000,0.230000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;];
%                     polyoffdefined = 0;
%                     cellstocutout = [];
                    
%                      % testing B from sweep 3 try 3
%                     final_matrix = [640,620,610,640,640,640,640,640,640,640,640,640,640,640;0.830000000000000,0.770000000000000,0.670000000000000,0.490000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000,0.450000000000000;0.230000000000000,0.410000000000000,0.650000000000000,0.920000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000,0.970000000000000;0.510000000000000,0.430000000000000,0.340000000000000,0.270000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000,0.260000000000000;]
%                     polyoffdefined = 0;
%                     cellstocutout = [];
% 
%                     % A with parallel (95%)
%                     final_matrix = [4110,4110,490,490,505,520,535,565,595,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685,685;0.950000000000000,0.950000000000000,0.900000000000000,0.850000000000000,0.800000000000000,0.1100000000000000,0.700000000000000,0.650000000000000,0.550000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000,0.400000000000000;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0.200000000000000,0.200000000000000,0.200000000000000,0.200000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000;];
%                     polyoffdefined = 1;
%                     obj.gratingDims.TL.polyoffset = .1470;
% 
%                     a_final = final_matrix(1,:);
%                     fill_final = final_matrix(2,:);
%                     r_final = final_matrix(3,:);
%                     off_final = final_matrix(4,:);
% %                     a_final(cellstocutout) = [];
% %                     fill_final(cellstocutout) = [];
% %                     r_final(cellstocutout) = [];
% %                     off_final(cellstocutout) = [];
% % %                     
% %                     obj.gratingDims.TL.GPbar = 10*round((a_final.*fill_final)/10)./1000;
% %                     a_round = 10*round(a_final/10)./1000;
% %                     obj.gratingDims.TL.GPgap = a_round-obj.gratingDims.TL.GPbar;
% %                     obj.gratingDims.TL.offset = 10*round(a_final.*off_final/10)./1000;
% %                     obj.gratingDims.TL.barpoly  = 10*round(r_final.*a_final.*fill_final/10)./1000;
% %                     
%                     obj.gratingDims.TL.GPgap = round(a_final.*(1-fill_final))./1000;
%                     obj.gratingDims.TL.GPbar = round(a_final.*fill_final)./1000;
%                     obj.gratingDims.TL.offset = round(a_final.*off_final)./1000;
%                     obj.gratingDims.TL.barpoly  = round(r_final.*a_final.*fill_final)./1000;
%                     if polyoffdefined ~= 1
%                         obj.gratingDims.TL.polyoffset = obj.gratingDims.TL.offset(1);
%                     end
% % %                     
%                     obj.gratingDims.TL.barpoly(1) = .040;
%                     obj.gratingDims.TL.barpoly(3) = .076;
%                     obj.gratingDims.TL.GPgap(1) = .1;
%                     obj.gratingDims.TL.GPgap(2) = .1;
                    
                    
%                     obj.gratingDims.TL.barpoly(1) = .04;
%                     obj.gratingDims.TL.GPgap(1) = .100;
%                     obj.gratingDims.TL.GPgap(2) = .100;

%                     % changes for best A ( BEST!!! )
%                     obj.gratingDims.TL.barpoly(3) = obj.gratingDims.TL.barpoly(3) - .002;
%                     obj.gratingDims.TL.offset(3) = obj.gratingDims.TL.offset(3) - .001;
%                     obj.gratingDims.TL.barpoly(2) = obj.gratingDims.TL.barpoly(2) + .001;

                    
%                     % final A layout
%                     layout = [0.505000000000000,0.404000000000000,0.404000000000000,0.126; 0.520000000000000,0.390000000000000,0.391000000000000,0.130; 0.536000000000000,0.3110000000000000,0.373000000000000,0.133; 0.565000000000000,0.367000000000000,0.367000000000000,0.141; 0.595000000000000,0.327000000000000,0.327000000000000,0.149; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; 0.685000000000000,0.200000000000000,0.200000000000000,0.110; ]

%                     % testing layout for gaussian shape
%                     layout = [.6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1; .6, .3, .3, .1;];
%                     obj.gratingDims.TL.polyoffset = .14;

%                     % CAPSTONE TESTING
%                      layout = [0.161000000000000,0.3110000000000000,0.3110000000000000,0.134000000000000;0.198000000000000,0.367000000000000,0.367000000000000,0.141000000000000;0.268000000000000,0.327000000000000,0.327000000000000,0.149000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;]
%                      obj.gratingDims.TL.polyoffset = 0.095
%                      layout(:,1) = layout(:,1) + layout(:,2);
%                      
%                      % CAPSTONE TESTING
%                      layout = [0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;0.250000000000000,0.200000000000000,0.200000000000000,0.110000000000000;];
%                      obj.gratingDims.TL.polyoffset = 0.095
%                      layout(:,1) = layout(:,1) + layout(:,2);
 
%                      % uniform
%                      layout = [0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;0.250,0.200,0.200,0.110;];
%                      obj.gratingDims.TL.polyoffset = 0.095
%                      layout(:,1) = layout(:,1) + layout(:,2);


                %%%%%%%%%%%%%%%% FINAL LAYOUTS V3 START %%%%%%%

%                 % final layout A (95% not DRC)
%                      layout = [0.475000000000000,0.451000000000000,0.451000000000000,0.0950; 0.475000000000000,0.451000000000000,0.451000000000000,0.0950; 0.490000000000000,0.441000000000000,0.441000000000000,0.0980; 0.491000000000000,0.417000000000000,0.417000000000000,0.0980; 0.505000000000000,0.404000000000000,0.404000000000000,0.126; 0.520000000000000,0.390000000000000,0.390000000000000,0.130; 0.536000000000000,0.375000000000000,0.375000000000000,0.134; 0.565000000000000,0.367000000000000,0.367000000000000,0.141; 0.595000000000000,0.327000000000000,0.327000000000000,0.149; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; 0.685000000000000,0.274000000000000,0.274000000000000,0.171; ]
%                      obj.gratingDims.TL.polyoffset = 0.095
                
                    % final layout A
                        layout = [0.510, 0.405, 0.405, 0.122;     0.525, 0.390, 0.390, 0.135;     0.540, 0.375, 0.375, 0.135;     0.570, 0.360, 0.360, 0.135;     0.600, 0.330, 0.330, 0.150;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     0.690, 0.270, 0.270, 0.165;     ]
                     obj.gratingDims.TL.polyoffset = 0.140
                    
%                     % final layout B
%                    layout = [0.649, 0.549, 0.0450, 0.360;     0.645, 0.540, 0.0750, 0.345;     0.634, 0.533, 0.0750, 0.345;     0.585, 0.405, 0.315, 0.180;     0.570, 0.360, 0.330, 0.150;     0.630, 0.300, 0.315, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     0.675, 0.255, 0.285, 0.150;     ]
%                     obj.gratingDims.TL.polyoffset = 1



                   %layout = round(layout*1000/10)/(1000/10);

                    obj.gratingDims.TL.GPgap = (layout(:,1) - layout(:,2)).';
                    obj.gratingDims.TL.GPbar = layout(:,2).';
                    obj.gratingDims.TL.offset = layout(:,4).';
                    obj.gratingDims.TL.barpoly = layout(:,3).';

                    %%%%%%%%%%%%%%%% FINAL LAYOUTS V3 END %%%%%%%

%                     

%                     % Changes for testing B sweep 2 try 2
%                     obj.gratingDims.TL.barpoly(1) = .04;
%                     obj.gratingDims.TL.GPgap(1) = .100;
%                     obj.gratingDims.TL.GPgap(2) = .100;
%                     obj.gratingDims.TL.offset(3) = .254 - .051;
%                     obj.gratingDims.TL.offset(4) = obj.gratingDims.TL.offset(4) + .079+.006;
%                     obj.gratingDims.TL.offset(5) = obj.gratingDims.TL.offset(5) + .079-.067;

                    
%                     % JUST KEEP TO CHANGE FROM OFFSET 2 to OFFSET1
%                     for j = 1:length(off_final)
%                         if off_final(j) < .5
%                             obj.gratingDims.TL.offset(j) = obj.gratingDims.TL.GPbar(j) - obj.gratingDims.TL.barpoly(j) - obj.gratingDims.TL.offset(j);
%                         else
%                             obj.gratingDims.TL.offset(j) = 2*obj.gratingDims.TL.GPbar(j) + obj.gratingDims.TL.GPgap(j) - obj.gratingDims.TL.barpoly(j) - obj.gratingDims.TL.offset(j);
%                         end
%                     end
                    
%                     obj.gratingDims.TL.GPgap =
%                     obj.gratingDims.TL.GPbar =
%                     obj.gratingDims.TL.offset =
%                     obj.gratingDims.TL.barpoly  = 
                    
                    
%%%%%%%%%%%%%%% CHANGES ENDING HERE!!!!!!!!
                  
                    %imperfections : added here by RK
                    obj.gratingDims.TL.widthSpacer = round(0.0/dx)*dx;
                    obj.gratingDims.TL.pSiImpW = round(0.0/dx)*dx;
                    obj.gratingDims.TL.pSiImpH = round(0.0/dx)*dx;
                    obj.gratingDims.TL.sinImpW = round(0.0/dx)*dx; %dip near body Si
                    obj.gratingDims.TL.sinImpH = round(0.0/dx)*dx; %dip near body Si
                    obj.gratingDims.TL.sin2ImpThickness = round(0.0/dx)*dx; %break up SiN liner into 2 materials
                    obj.gratingDims.TL.sin2_index_1550 = 2.3;
                otherwise
                    error('No matching design. Please see help file for valid design inputs')
            end
            
        end
        
        function obj = varLoops(obj)
            if isfield(obj.gratingDesign,'varLoops')
                if isfield(obj.gratingDesign.varLoops,'TL')
                    varFNTL = fieldnames(obj.gratingDesign.varLoopsBase.TL);
                    switch obj.gratingDesign.varLoops.varType
                        case 'percent'
                            for ii = 1:length(varFNTL)
                                eval(['obj.gratingDims.TL.' varFNTL{ii} '=obj.gratingDims.TL.' varFNTL{ii} '*(1+obj.gratingDesign.varLoops.TL.' varFNTL{ii} '_ITER)'])
                            end
                        case 'absolute'
                            for ii = 1:length(varFNTL)
                                eval(['obj.gratingDims.TL.' varFNTL{ii} '= obj.gratingDesign.varLoops.TL.' varFNTL{ii} '_ITER'])
                            end
                        otherwise
                            error('Specify type of variation. percent or absolute')
                    end
                end
                
            end
        end
        
    end
    
end

