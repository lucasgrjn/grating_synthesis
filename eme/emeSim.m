%%
% Matlab class 'emeSim'
%
% Authors: Cale Gentry
%
% All spatial units in microns


%%
classdef emeSim
    
    properties
        constants;                                                              %Constants (speed of light, impedance of free space)
        domain;                                                                 %Properties of the domain (size, discretization, etc...)
        structures;                                                             %Structures drawn in the domain (layers, rectangles, etc...)
        diel;                                                                   %The dielectric profile itself (refractive index, not relative permittivity)
        layerProperties;
        modes;
        scatterProperties;
        fullFields;
        sources;
        fiberCoup;
        obsPlanes;
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%              Set Up emeSim Class                 %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = emeSim(varargin)
            
            obj.constants.c = 299792458;                                        %SI units [m/s]
            obj.constants.mu0 = 4e-7*pi;                                        %SI units [m kg/(s^2 A^2)]
            obj.constants.eps0 = 1/(obj.constants.c^2*obj.constants.mu0);       %SI units [s^4 A^2/(m^3 kg)]
            
            domainFields = {'table', 'domain', ...                              %List of domain fields allowed to be passed to emeSim()
                'discretization','pml','wavelengthSpectrum','debug', ...
                'polarization', 'backgroundIndex'};
            for k = 1:2:size(varargin,2)                                        %Assigns the value after each field to be added to obj.domain. structure
                varName = varargin{k};  
                switch(varName)
                    case domainFields
                        eval(['obj.domain.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid' ...
                            ' parameter "' varargin{k} ...
                            '".  Check fmm instructions.']);
                end
            end
            
            % Set defaults for unset values or throw errors
            if ~isfield(obj.domain,'discretization')                             
                error('You must specify the discretization');
            end
            if ~isfield(obj.domain,'domain')
                error('You must specify the dimensions of the domain');
            else
                xf = obj.domain.domain(1); zf = obj.domain.domain(2);           %Assign the domain values to obj.domain.domain vector
                dx = obj.domain.discretization(1); 
                dz = obj.domain.discretization(2);
                if round(xf/dx) ~= round(xf/dx*10^6)/10^6                       %Checks to see if domain is integer of specified discretization 
                    fprintf(['\n Warning: x domain is not an integer number'...
                        'of discretization points \n']); 
                end
                if round(zf/dz) ~= round(zf/dz*10^6)/10^6
                    fprintf(['\n Warning: z domain is not an integer number'...
                        'of discretization points \n']);
                end
                if ~isfield(obj.domain,'backgroundIndex')
                    obj.domain.backgroundIndex = 1;
                    fprintf(['\n Warning: no backgroundIndex input.'...
                        ' Defaulting to air (n = 1) \n']);
                end
                obj.diel = ones(round(xf/dx), round(zf/dz))*...
                    obj.domain.backgroundIndex;
            end
            if ~isfield(obj.domain,'pml')                                       %PML information is required to be assigned
                error('You must specify the thicknesses of the pmls');
            end
            if ~isfield(obj.domain,'wavelengthSpectrum')                        %Must specify wavelength
                error('wavelengthSpectrum was not set.');
            end
            if length(obj.domain.pml) == 1
                fprintf(['\nOnly thickness of PML set. Setting to quadratic '...
                    'with xim = 20.\n']);
                obj.domain.pml(2) = 2;
                obj.domain.pml(3) = 20;
            elseif length(obj.domain.pml) == 2
                fprintf(['\nDistance of pml into imaginary plane not set '...
                    'setting xim = 20.\n']);
                obj.domain.pml(3) = 20;
            end
            if obj.domain.pml(2) < 1
                error('\nInvalid pml order (must be greater than 0\n');
            end
            if length(obj.domain.wavelengthSpectrum) == 1
               lambda = obj.domain.wavelengthSpectrum;
               fprintf(['\n Single wavelength of ' num2str(lambda) ' \mum \n']);
               obj.domain.wavelengthSpectrum = [lambda lambda 1];
            end
            if ~isfield(obj.domain,'debug')                                     %This might not end up being used
                obj.domain.debug = 'no';
                fprintf('\nDebug was not specified.  Setting to no.\n');
            end
            if ~isfield(obj.domain,'polarization')                              %Must assign polarization
                error(['polarization was not set. You must specify the'...
                    ' polarization (''TE'' or ''TM'') for a 2D simulation.\n']);
            end
           
            % update z for consistency
            % NOTE: dz is recalculated here and may differ from what the
            % user intended
            % domain expects z to start at dz and end at zf
            obj.domain.z = linspace(0, obj.domain.domain(2), size(obj.diel, 2)+1 );
            obj.domain.z = obj.domain.z(2:end);
            obj.domain.discretization(2) = obj.domain.z(2) - obj.domain.z(1);
            
            %{
            %future
            if ~isfield(obj.domain,'colorScheme')
                obj.domain.colorScheme = 'redbluedark';
                fprintf('colorScheme was not specified. Setting to redWhiteBlue.\n');
            end
            
                
            if ~isfield(obj.domain,'namedd')
                obj.domain.namedd = 'scatteringMatrix_dd';
                fprintf('namedd was not set. Scattering matrix will be saved in the text file ''scatteringMatrix_dd.txt''.\n');
            end
            if ~isfield(obj.domain,'nameDir')
                obj.domain.nameDir = 'fdtdSimulationData';
                fprintf('nameDir was not set. The data files will be transferred to the directory ''fdtdSimulationData''.\n');
            end
            %}
        end                           %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%        Add Planar Layer(s) To The Domain         %-------%%%
        %%%----------------------------------------------------------%-------%%%
        function obj = addLayers(obj,varargin)
            
            fprintf('\nAdding layers. \n');
            preNumStruc = length(obj.structures);                               %Get the number of structures already defined
            layerFields = {'numLayers','centerZ','widthZ','minX','heightX',...
                'indices'};
            
            for k = 1:2:size(varargin,2)                                        %Assign the values to the layerFields
                varName = varargin{k}; 
                switch(varName)
                    case layerFields
                        eval(['temp.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter'...
                            ' "' varargin{k} '".  Check fmmSim instructions.']);
                end
            end
            
            if ~isfield(temp,'numLayers')                                       %If number of layers isn't assigned, throw an error
                error(['You must specify a number of layers when calling'...
                    ' addLayers. Check fmmSim instructions']);
            end
            
            if ~isfield(temp,'centerZ')                                         %If the centerZ of the layer isn't defined, default to the middle of the computational domain
                fprintf(['\nCenter z-coordinates of layers not set.  Setting'...
                    ' to default value of the middle of the'...
                    ' simulation space.\n']);
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).centerZ = ...
                        obj.domain.domain(2)/2;
                end
            elseif ( length(temp.centerZ) == temp.numLayers )
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).centerZ = temp.centerZ(ii);
                end
            else
                error('centerZ must have length numLayers.');                   %Input vector of centerZ's and numLayers must be the same
            end
            
            if ~isfield(temp,'widthZ')                                          %If the width of the layer isn't defined default to the entire width of the computational domain
                fprintf(['Widths in z of layers not set.  Setting length of'...
                    ' computational space \n']);
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).widthZ = ...
                        obj.domain.domain(2);
                end
            elseif ( length(temp.heightX) == temp.numLayers )
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).widthZ = temp.widthZ(ii);
                end
            else
                error('widthZ must have length numLayers.');                    %Can't have widthZ a different size than numLayers (Note: This means you can't have one layer default while defining another)
            end
            
            if ~isfield(temp,'minX')                                            %If the bottom of the layer isn't defined, throw error
                error(['You must specify the minimum x values of all the'...
                    'layers to be added to the structure. Check fmmSim'...
                    'instructions.']);
            elseif ( length(temp.minX) == temp.numLayers )
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).minX = temp.minX(ii);
                end
            else
                error('minX must have length numLayers.');
            end
            
            if ~isfield(temp,'heightX')                                         %If the height of the layer isn't defined, throw error
                error('Height of x layers is not set.'); 
            elseif ( length(temp.heightX) == temp.numLayers )
                for ii = 1:temp.numLayers
                    obj.structures(preNumStruc+ii).heightX = temp.heightX(ii);
                end
            else
                error('heightX must have length numLayers.');
            end
            
            %Get size and discretization of domain
            xf = obj.domain.domain(1); zf = obj.domain.domain(2); 
            dx = obj.domain.discretization(1); 
            dz = obj.domain.discretization(2);
            
            %Assign equivalent anchor point, structureDimensions, and index
            for ii = 1:temp.numLayers                                           
                obj.structures(preNumStruc+ii).anchorPoint = ...                %Assign all the inputs to obj.structures(structure #)
                    [obj.structures(preNumStruc+ii).minX ...
                    (obj.structures(preNumStruc+ii).centerZ-obj.structures(...
                    preNumStruc+ii).widthZ/2.)];
                obj.structures(preNumStruc+ii).structureDimensions = ...
                    [obj.structures(preNumStruc+ii).heightX ...
                    obj.structures(preNumStruc+ii).widthZ];
                obj.structures(preNumStruc+ii).index = temp.indices(ii);
                ancx = obj.structures(preNumStruc+ii).minX; ancz = ...
                    (obj.structures(preNumStruc+ii).centerZ-obj.structures(...
                    preNumStruc+ii).widthZ/2.);
                height = obj.structures(preNumStruc+ii).heightX; lengthz = ...
                    obj.structures(preNumStruc+ii).widthZ; 
                writestruc = 'yes';                                             %Initialize writestruc to yes. This will get turned to no if you wanted to make a structure out of the computational domain.
                
                %X parameters (vertical)
                if round(ancx/dx)~= round(ancx/dx*10^6)/10^6                    %Check to see if it's assigned on the grid  (the 10^6 makes it so there isn't numerical precision problems)
                   fprintf(['\n Warning: Bottom of layer Structure ' num2str(...
                       preNumStruc + ii) ' is not a integer of dx \n']); 
                end
                %If layers starts below computational domain
                if round(ancx/dx) < 1   
                    if round((ancx+height)/dx) < 1                              %If it ends under the computatonal domain also
                        fprintf(['\n Warning: Bottom of layer Structure '...
                            num2str(preNumStruc + ii)...
                            ' is below the domain \n']); 
                        writestruc = 'no';
                    else startx = 1;                                            %Start writing on the first discretation point
                    end
                else startx =  round(ancx/dx)+1;                                %ie if the layer starts at 1 um and the discretization is .1 um then that index gets put on the 11th pt (just like 2DFDTD)
                end
                
                if round((ancx+height)/dx) > round(xf/dx)                       %If the top of the structure is definied out of the domain, truncate it.
                    if round((ancx)/dx) > round(xf/dx)                             %Unless all of it is out of the domain then don't write it and give a warning
                        fprintf(['\n Warning: Layer Structure' num2str(...
                            preNumStruc + ii) 'is above the domain \n']); 
                        writestruc = 'no';
                    else endx = round(xf/dx); 
                    end
                else endx = round((ancx+height)/dx);
                end
                
                %Z parameters (horizontal... propagation direction)
                if round(ancz/dz)~= round(ancz/dz*10^6)/10^6                    %If the front of the layer isn't on the grid, give warning
                   fprintf(['\n Warning: Layer Structure ' num2str(...
                       preNumStruc + ii) ' width is not a integer of dz \n']); 
                end
                if round(ancz/dz) < 1                                           %If the left side of the layer is to the left of the computational domain
                    if round((ancz+lengthz)/dz) < 1                                 %If the right side of the layer is also left of the computational domain, don't write the layer
                        fprintf(['\n Warning: Layer Structure ' num2str(...
                            preNumStruc + ii) ' is left of the domain \n']); 
                        writestruc = 'no';
                    else startz = 1;
                    end
                else startz =  round(ancz/dz);   
                end
                if round((ancz+lengthz)/dz) > round(zf/dz)                      %If the right side of the layer is to the right of the computational domain
                    if round((ancz)/dz) > round(zf/dz)                               %If the left side of the layer is also to the right of the computational domain, don't write the layer
                        fprintf(['\n Warning: Structure' num2str(...
                            preNumStruc + ii) 'is right of the domain \n']); 
                        writestruc = 'no';
                    else endz = round(zf/dz);
                    end
                else endz = round((ancz+lengthz)/dz);
                end
                
                %Actually draw structure
                if ~strcmpi(writestruc,'no')
                    obj.diel(startx:endx, startz:endz) = obj.structures(...
                        preNumStruc + ii).index;
                end
            end
        end                    %-------%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%          Add Shape(s) To The Domain              %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = addStructure(obj,varargin)
            %
            % Example:
            % Q = Q.addStructure('Rectangle', ...
            %     'anchorPoint',[(hAir+hBotOx) start], ...  %Anchor POint is at the [bottom left]
            %     'numRepeat', numPeriods, ...
            %     'repeatPeriod', period, ...
            %     'repeatVec', [0 0 1], ...                        %Repeat in the z-direction (almost always this way)
            %     'structureDimensions',[hSi gap], ...       %[height 
            %     'index', nSiO2);
            
            
            objectTypes = {'Rectangle','Triangle','RingSector',...              %Currently only rectangle works
                'EllipticSector','ObjectFromFile2D'};
            varName = varargin{1};
            fprintf(['\nAdding a ',varName,' \n'])
            switch (varName)
                case objectTypes
                    idx = strcmpi(varName,objectTypes);
                    idx = find( idx == 1 );                                     
                otherwise
                    error(['Trying to add undefined structure type "'...
                        varargin{1} '".  Check EME instructions by typing '...
                        'help emeSim''.']);
            end
            
            numStruc = length(obj.structures) + 1;
            if idx == 1%rectangle                                               %If rectangle create structure fields
                structureFields = {'anchorPoint','structureDimensions',...
                    'orientationVecs','index','numRepeat','repeatPeriod',...
                    'repeatVec'};
            elseif idx > 1                                                      %If not... throw error
                error(['Sorry, the geometric structure "' objectTypes{idx}...
                    '" is not supported at this time. Contact Cale to make'...
                    ' it happen.'])
            end
            for k = 2:2:size(varargin,2)                                        %Assign values to structuresFields
                varName = varargin{k};
                switch(varName)
                    case structureFields
                        eval(['obj.structures(numStruc).'...
                            varName ' = varargin{k+1};']);
                    otherwise                                                   %Throw error if they assigned a structureField that doesnt exist
                        error(['Trying to assign value to invalid'...
                            ' parameter "' varargin{k} '".  Check EME'...
                            ' instructions by typing ''help emeSim''.']);
                end
            end
            if (~isfield(obj.structures(numStruc),'numRepeat') ||...            %If numRepeat doesn't exist OR it has no definied valued, set to one
                    isempty(obj.structures(numStruc).numRepeat) )
                fprintf('numRepeat not set.  Setting to 1.\n');
                obj.structures(numStruc).numRepeat = 1;
            elseif obj.structures(numStruc).numRepeat > 1                       %If numRepeat is set to something greater than 1
                if or(~isfield(obj.structures(numStruc),'repeatPeriod'), ...    %  if it repeatPeriod isn't set, throw an error
                        isempty(obj.structures(numStruc).repeatPeriod) )
                    error(['For a repeating structure you must specifiy the'...
                        ' period of repetition.']);
                elseif or( ~isfield(obj.structures(numStruc),'repeatVec'), ...  %  or if repeatVec isn't set, throw an error
                        isempty(obj.structures(numStruc).repeatVec) )
                    error(['For a repeating structure you must specifiy the'...
                        ' vector along which repetition occurs.']);
                end
            end
            if or( ~isfield(obj.structures(numStruc),'anchorPoint'),...         %If anchorPoint isn't set, throw error
                    isempty(obj.structures(numStruc).anchorPoint) )
                error(['You must specify the anchorPoint for structure '...
                    num2str(numStruc)]);
            end
            if length(obj.structures(numStruc).anchorPoint) ~= 2                %If anchor point doesnt contain 2 numbers, throw error
                error(['Your anchor coordinate point is not valid. Must be'...
                    ' of form [x-coordinate z-coordinate] for structure ' ...
                    num2str(numStruc)]);
            end
            if or( ~isfield(obj.structures(numStruc),'structureDimensions'),... %If structure dimensions aren't set, throw error
                    isempty(obj.structures(numStruc).structureDimensions) )
                error(['You must specify the structureDimensions for'...
                    ' structure ' num2str(numStruc)]);
            end
            if length(obj.structures(numStruc).structureDimensions) ~= 2        %If structure dimensions doesnt contain 2 numbers, throw error
                error(['Your structure dimensions is not valid. Must be of'...
                    ' form [x-length z-length] for structure '...
                    num2str(numStruc)]);
            end
            if ( ~isfield(obj.structures(numStruc),'index') || ...              %If the refractive index of the structure isn't set, throw error
                    isempty(obj.structures(numStruc).index) )
                error(['You must specify the index of refraction for'...
                    ' structure ' num2str(numStruc)]);
            end         
            
            %Get domain size and discretization
            xf = obj.domain.domain(1); zf = obj.domain.domain(2); 
            dx = obj.domain.discretization(1); 
            dz = obj.domain.discretization(2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %--------------------------RECTANGLE-------------------------------%
            if idx == 1 
                height = obj.structures(numStruc).structureDimensions(1);
                lengthz = obj.structures(numStruc).structureDimensions(2);
                for ii = 1:obj.structures(numStruc).numRepeat                   %Assign defined structure properties to obj.strcutures(structure #) 
                    if obj.structures(numStruc).numRepeat == 1
                        ancx = obj.structures(numStruc).anchorPoint(1); 
                        ancz = obj.structures(numStruc).anchorPoint(2);
                    elseif isequal(obj.structures(numStruc).repeatVec, [0 0 1]) %If the period is defined in the z direction
                        ancx = obj.structures(numStruc).anchorPoint(1); 
                        if round(obj.structures(numStruc).repeatPeriod/dz)...   %If the period isn't an integer of the discretization, give warning
                                ~= round((obj.structures(numStruc).repeatPeriod...
                                /dz)*1e12)*1e-12
                            fprintf(['\n' num2str(obj.structures(...
                                numStruc).repeatPeriod) ' Period not integer'...
                                ' number of discretizations \n']);
                        end
                        ancz = obj.structures(numStruc).anchorPoint(2)+...      %Shifts the anchor point depending on what period it is writing
                            round(obj.structures(numStruc).repeatPeriod/dz)*...
                            dz*(ii-1);
                    elseif isequal(obj.structures(numStruc).repeatVec, [1 0 0]) %If the period is defined in the x direction
                        if round(obj.structures(numStruc).repeatPeriod/dx)~=...
                                (obj.structures(numStruc).repeatPeriod/dx)
                            fprintf(['\n' num2str(obj.structures(numStruc...    %Throw warning if the period isn't an integer of discretization.
                                ).repeatPeriod) ' Period not integer number'...
                                'of discretizations \n']);
                        end
                        ancx = obj.structures(numStruc).anchorPoint(1)+...      %Shifts the anchor point depending on what period it is writing
                            round(obj.structures(numStruc).repeatPeriod/dx)*...
                            dx*(ii-1); 
                        ancz = obj.structures(numStruc).anchorPoint(2);
                    else error(['repeatVec = ', num2str(obj.structures(...      %Throw an error if the direction of the repeatVec isn't x or z
                            numStruc).repeatVec), ' is not valid. It must be'...
                            ' [1 0 0] or [0 0 1]']);
                    end
                    if ii == 1
                        if round(ancx/dx) ~= round(ancx/dx*10^6)/10^6           %Check to see if the bottom of the rectangle is on the grid, else give warning
                            fprintf(['\n Warning: x anchor point is not an'...
                                ' integer number of discretization points \n']); 
                        end
                        if round(ancz/dz) ~= round(ancz/dz*10^6)/10^6           %Check to see if the left side of the rectangle is on the grid, else give warning
                            fprintf(['\n Warning: z anchor point is not an'...
                                ' integer number of discretization points \n']);
                        end
                        if round(height/dx) ~= round(height/dx*10^6)/10^6       %Check to see if the height (x dimension) of the rectangle is on the grid, else give warning
                            fprintf(['\n Warning: x dimension is not an'...
                                ' integer number of discretization points \n']); 
                        end         
                        if round(lengthz/dz) ~= round(lengthz/dz*10^6)/10^6     %Check to see if the length (z dimension) of the rectangle is on the grid, else give warning
                            fprintf(['\n Warning: z dimension is not an'...
                                'integer number of discretization points \n']); 
                        end
                    end
                    writestruc = 'yes'; %initialize writestruc
                    
                    %X parameters (vertical)
                    if round(ancx/dx) < 1                                       %If the bottom of the rectangle is below the computational domain
                        if round((ancx+height)/dx) < 1                              %If the top of the rectangle is also below the computational domain, give warning and don't write structure
                            fprintf(['\n Warning: Structure '...
                                num2str(numStruc) ' period ' num2str(ii)...
                                'is below the domain \n \n']); 
                            writestruc = 'no';
                        else startx = 1;
                        end
                    else startx =  round(ancx/dx)+1;                            
                    end 
                    if round((ancx+height)/dx) > round(xf/dx)                   %If the top of the rectangle is above the computational domain           
                        if round((ancx)/dx) > round(xf/dx)                           %If the bottom of the rectangle is also above the computational domain, give warning and don't write structure
                            fprintf(['\n Warning: Structure '...
                                num2str(numStruc) ' period ' num2str(ii)...
                                ' is above the domain \n \n']); 
                            writestruc = 'no';
                        else endx = round(xf/dx);
                        end
                    else endx = round((ancx+height)/dx);
                    end
                    
                    %Z parameters (horizontal... propagation direction)
                    if round(ancz/dz) < 1                                       %If the left side of the rectangle is left of the computaitonal domain
                        if round((ancz+lengthz)/dz) < 1                             %If the right side of the rectangle is also to the left of the computational domain, give warning and don't write structure
                            fprintf(['\n Warning: Structure '...
                                num2str(numStruc)  ' period ' num2str(ii)...
                                ' is left of the domain \n \n']); 
                            writestruc = 'no';
                        else startz = 1;
                        end
                    else startz =  round(ancz/dz)+1;   
                    end
                    if round((ancz+lengthz)/dz) > round(zf/dz)                  %If the right side of the rectangle is to the right of the computational domain
                        if round((ancz)/dz) > round(zf/dz)                          %If the left side of the rectangle is also to the right of the computational domain, give warning and don't write structure
                            fprintf(['\n Warning: Structure '...
                                num2str(numStruc)  ' period ' num2str(ii)...
                                ' is right of the domain \n \n']); 
                            writestruc = 'no';
                        else endz = round(zf/dz);
                        end
                    else endz = round((ancz+lengthz)/dz);    
                    end
                    
                    %Actually draw structure
                    if ~strcmpi(writestruc,'no')
                        obj.diel(startx:endx, startz:endz) = obj.structures(...
                            numStruc).index;
                    end
                end %each period of rectangle
            end %if for rectangle           
        end                 %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%           Convert Dielectric Matrix              %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = convertDiel(obj)
            fprintf('\nCompressing dielectric information.\n');
            
            %z dimension (propagation direction)
            zloc = find(1-(all(diff(obj.diel,1,2)==0)))+1;                      %List of z indeces where interfaces are. These indeces are of the first new layer. (-1 would last of previous layer)
            listOfLayers = obj.diel(:,[1 zloc]);                                %All the layer index distributions
            
            
            zOfLayers = [zloc size(obj.diel,2)]*...                             %The z position of each layer
                obj.domain.discretization(2);
            lengthOfLayers = [zOfLayers(1,1) diff(zOfLayers)];                  %This is a list of all the z thicknesses, except the first and last that will be one discretization off
                lengthOfLayers(1,1) = lengthOfLayers(1,1)-...
                    obj.domain.discretization(2);                                   %Shift first one
                lengthOfLayers(1,end) = lengthOfLayers(1,end)+...
                    obj.domain.discretization(2);                                   %Shift last one
            
            [uniqueCrossSections, ia, crossSectionNum] = ...                    %Pick out the unique layers to find the modes of and store which ones they are 
                unique(listOfLayers.', 'rows', 'stable');                       %   'stable' outputs them in the order they are found. Not uniqueCrossSections are in the rows here
            
            obj.layerProperties.numOfLayers = size(listOfLayers,2);             %Number of layers
            obj.layerProperties.uniqueCrossSections = uniqueCrossSections.';    %Assign them to the obj.layerProperties. We want uniqueCrossSections to be columns
            obj.layerProperties.crossSectionNum = crossSectionNum;              %This gives the corresponding unique layer for a given actual layer (ie crossSectionNum(layerNum) = whichUniqueCrossSection)
            obj.layerProperties.numOfUniqueCrossSections = length(ia);          %Number of unique cross sections
            obj.layerProperties.lengthOfLayers = lengthOfLayers;                %Length of each layer
            obj.layerProperties.numOfInterfaces = ...                           %Number of interfaces
                obj.layerProperties.numOfLayers-1;
            
            interfaces = [obj.layerProperties.crossSectionNum(1:end-1,1)...     %List of all the interfaces where left row is left layer and right row is right layer
                obj.layerProperties.crossSectionNum(2:end,1)];
            [obj.layerProperties.uniqueInterfaces, ia2, ...                     %Assign the uniqueInterfaces where the left row corresponds to the left layer defined by which unique layer it is and likewise with the right 
                obj.layerProperties.interfaceNum] = ...
                unique(interfaces, 'rows', 'stable');                           %   'stable' outputs them in the order they are found
            obj.layerProperties.numOfUniqueInterfaces = length(ia2);            %Number of unique interfaces
            
            
            xf = obj.domain.domain(1); zf = obj.domain.domain(2); 
            dx = obj.domain.discretization(1); 
            dz = obj.domain.discretization(2);
            [xSize, zSize] = size(obj.diel);
            if round(xf/dx*1e6)*1e-6 ~= xSize
                fprintf(['\n Warning: xf parameter and dx give a different ' ...
                    'size computational domain than Q.diel.  This can cause errors in addObsPlane \n']); 
            end
            
            if round(zf/dz*1e6)*1e-6 ~= zSize
                fprintf(['\n Warning: zf parameter and dz give a different ' ...
                    'size computational domain than Q.diel.  This can cause errors in addObsPlane \n']); 
            end
            
                
            
        end                           %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%            Add Source Information                %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = addSource(obj,varargin)
            
            fprintf('\nAdding a source.\n');

            sourceFields = {'modeNumber','plotSource'};                         %The modenumber (in order given by modeslv) and whether or not to plot the source
            for k = 1:2:size(varargin,2)                                        %Assign the values passed 
                varName = varargin{k};
                switch(varName)
                    case sourceFields
                        eval(['obj.sources.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid'...
                            ' parameter "' varargin{k} '".  Check EME'...
                            ' instructions by typing ''help emeSim''.']);
                end
            end
            if ~isfield(obj.sources,'modeNumber')                               %If modeNumber isn't set, set to 1, which is usually the fundamental mode
                obj.sources.modeNumber = 1;
                fprintf(['modeNumber was not set.  Setting to default value'...
                    '1.\n']);
            end
            if ~isfield(obj.sources,'plotSource')                               %If plotSource isn't set, set it to no and print that it did so.
                obj.sources.plotSource = 'no';
                fprintf('plotSource was not set. Setting to no.\n');
            end
            if strcmpi(obj.sources.plotSource, 'yes')
                source = zeros(length(obj.modes.fields{1,1}(:,1,1)),...
                    length(obj.domain.k0));
                for ii = 1:length(obj.domain.k0)
                    source(:,ii) = ...
                        obj.modes.fields{1,1}(:,obj.sources.modeNumber,ii);
                end
                
                if obj.domain.polarization == 0
                    [temp, xaxis] = meshgrid(1:length(obj.domain.k0),...
                    obj.domain.discretization(1):obj.domain.discretization(1)...
                    :obj.domain.domain(1));
                elseif obj.domain.polarization == 1
                    [temp, xaxis] = meshgrid(1:length(obj.domain.k0),...
                    obj.domain.discretization(1):obj.domain.discretization(1)...
                    :obj.domain.domain(1)-obj.domain.discretization(1));
                end

                figure; 
                plot(xaxis, real(source), 'LineWidth', 2); grid on; 
                set(gcf, 'Color', 'White'); 
                set(gca, 'FontSize', 14);
                xlabel('x (\mum)','FontSize',14);
                ylabel('Real Field (AU)','FontSize',14); 
                title('Source Modes', 'FontSize', 16);
                legend([num2str((obj.domain.wavelengthSpectrum(1):...
                    obj.domain.wavelengthSpectrum(3):...
                    obj.domain.wavelengthSpectrum(2))'*10^3)]);
            end
        end                    %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%          Find The Modes For Each Layer           %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = getModes(obj, n_modes)
            %Still need to take the defaults off the PML order and depth
            % Inputs:
            %   n_modes
            %       type: integer
            %       desc: optional input for # of modes to solve for.
            %               Default is solve for all the modes
            
            % get number of modes
            if nargin > 1
                OPTS.n_modes = n_modes;
            end
            
            obj.domain.k0 = 2*pi./(obj.domain.wavelengthSpectrum(1):...         %Create the free space wavenumber vector
                obj.domain.wavelengthSpectrum(3):...
                obj.domain.wavelengthSpectrum(2));
            pol = obj.domain.polarization;                                      %Get the polarization
            
            %Get PML information
            OPTS.PML.d = obj.domain.pml(1);                                        %Width of PML (microns)
            OPTS.PML.m = obj.domain.pml(2)-1;                                                     %Order of PML (1 means quadratic)
            OPTS.PML.xim = obj.domain.pml(3);                                                  %Length of PML into imaginary space (mirons)

            dx = obj.domain.discretization(1);                                  %Only need x discretization for this section

            obj.modes.fields = cell(...                                         %Places to store information about each layer
                obj.layerProperties.numOfUniqueCrossSections,1);                %The fields Ey for TE and Hy for TM
            B = cell(obj.layerProperties.numOfUniqueCrossSections,1);           %The material parameters (intermediate mus (TE) or eps (TM)) also includes PML info
            obj.modes.betas = cell(...
                obj.layerProperties.numOfUniqueCrossSections,1);                %The propagation constants
            fundbetas = cell(obj.layerProperties.numOfUniqueCrossSections,1);   %The propatation constant of the fundamental modes
            obj.modes.qs = cell(obj.layerProperties.numOfUniqueCrossSections,1);%The row/column of the fundamental mode in each layer (r from modeslv)
            n = cell(obj.layerProperties.numOfUniqueCrossSections,1);

            for j = 1:obj.layerProperties.numOfUniqueCrossSections              %Run modesolver in each unique layer
                for ii = 1:length(obj.domain.k0)
                [obj.modes.fields{j}(:,:,ii),  obj.modes.betas{j}(:,:,ii),...
                    B{j}(:,:,ii), fundbetas{j}(ii),obj.modes.qs{j}(ii), n{j}]...
                    = modeslv(obj.layerProperties.uniqueCrossSections(:,j), ...
                    dx, obj.domain.k0(ii), pol, OPTS);
                end 
            end
            
            obj.modes.B = B;                                                    %Store B 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Overlap matrices at each interface
            %|O12 O21|
            %|O23 O32|
            %|O34 O43|
            %|O45 O54|
            %| .   . |
            %| :   : |
            obj.layerProperties.overs = cell(...                                %Initialize Os
                obj.layerProperties.numOfUniqueInterfaces,2);
            for j = 1:obj.layerProperties.numOfUniqueInterfaces
                for ii = 1:length(obj.domain.k0)
                    obj.layerProperties.overs{j,1}(:,:,ii) = ...                %Overlap from left to right   FL.'*B*FR
                        .5*1/((obj.domain.k0(ii)*10^6*obj.constants.c))*...
                        ((obj.modes.fields{...
                        obj.layerProperties.uniqueInterfaces(j,1)}(:,:,ii)...
                        ).')*B{obj.layerProperties.uniqueInterfaces(j,2)}(:,:...
                        ,ii)*obj.modes.fields{...
                        obj.layerProperties.uniqueInterfaces(j,2)}(:,:,ii)*...
                        obj.modes.betas{obj.layerProperties.uniqueInterfaces(...
                        j,2)}(:,:,ii)*...
                        obj.domain.discretization(1)*10^-6;
                    obj.layerProperties.overs{j,2}(:,:,ii) = ...                %Overlap from right to left  FR.'*B*FL
                        .5*1/((obj.domain.k0(ii)*10^6*obj.constants.c))*...
                        ((obj.modes.fields{...
                        obj.layerProperties.uniqueInterfaces(j,2)}(:,:,ii)...
                        ).')*B{obj.layerProperties.uniqueInterfaces(j,1)}(:,:...
                        ,ii)*obj.modes.fields{...
                        obj.layerProperties.uniqueInterfaces(j,1)}(:,:,ii)*...
                        obj.modes.betas{obj.layerProperties.uniqueInterfaces(...
                        j,1)}(:,:,ii)*...
                        obj.domain.discretization(1)*10^-6; 
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Transfer matrices at each interface
            %|T12 T21|
            %|T23 T32|
            %|T34 T43|
            %|T45 T54|
            %| .   . |
            %| :   : |
            obj.layerProperties.Ts = cell(...                                   %Initialize Ts
                obj.layerProperties.numOfUniqueInterfaces,2); 
            for j = 1:obj.layerProperties.numOfUniqueInterfaces
                for ii = 1:length(obj.domain.k0)
                    obj.layerProperties.Ts{j,1}(:,:,ii) =...                    %T12 = 2*inv(O12+O21.')
                        2*inv(obj.layerProperties.overs{j,1}(:,:,ii)+...
                        (obj.layerProperties.overs{j,2}(:,:,ii)).');
                    obj.layerProperties.Ts{j,2}(:,:,ii) = ...                   %T21 = T12.'
                        (obj.layerProperties.Ts{j,1}(:,:,ii)).';
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Reflection matrices at each interface
            %|R12 R21|
            %|R23 R32|
            %|R34 R43|
            %|R45 R54|
            %| .   . |
            %| :   : |
            obj.layerProperties.Rs = cell(...                                   %Initialize Rs
                obj.layerProperties.numOfUniqueInterfaces,2);  
            for j = 1:obj.layerProperties.numOfUniqueInterfaces
                for ii = 1:length(obj.domain.k0)
                    obj.layerProperties.Rs{j,1}(:,:,ii) =...                    %R12 = 1/2(O21.'-O12)*T12
                        0.5*((obj.layerProperties.overs{j,2}(:,:,ii)).'-...
                        obj.layerProperties.overs{j,1}(:,:,ii))*...
                        obj.layerProperties.Ts{j,1}(:,:,ii);
                    obj.layerProperties.Rs{j,2}(:,:,ii) =...                    %R21 = 1/2(O12.'-O21)*T21
                        0.5*((obj.layerProperties.overs{j,1}(:,:,ii)).'-...
                        obj.layerProperties.overs{j,2}(:,:,ii))*...
                        obj.layerProperties.Ts{j,2}(:,:,ii);
                end
            end
            
        end                              %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%      Create Scatter Matrices and Fields          %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = scatter(obj)
            %Needs fix:Things might mess up for 3 or less layers
            
            %if any fields clear them
            if isfield(obj.scatterProperties, 'PowerIn')
                obj.scatterProperties = ...
                    rmfield(obj.scatterProperties,fieldnames(obj.scatterProperties));
            end
            
            %Create t and r matrices
            numodes = length(obj.modes.betas{1,1}(:,1,1));                      %Get the number of modes (also the number of pixels in the x direction)
            ts = cell(obj.layerProperties.numOfLayers, ...                      %Instantiate cells
                obj.layerProperties.numOfLayers);
            rs = cell(obj.layerProperties.numOfLayers, ...
                obj.layerProperties.numOfLayers);

            fprintf('Computing S matrices and propagating...\n');
            
            tic;    % DEBUG
            for jj = 1:obj.layerProperties.numOfInterfaces                       %For number of interfaces
                for ii = 1:length(obj.domain.k0)                                    %For each wavelength
                    
%                     % OLD VERSION, much more inefficient
%                     % transfer 
%                     ts{jj,jj+1}(:,:,ii) = obj.layerProperties.Ts{ obj.layerProperties.interfaceNum(jj), 1 }(:,:,ii) * ...      
%                         expm( -1i* obj.modes.betas{ obj.layerProperties.crossSectionNum(jj, 1) }(:,:,ii) * ...
%                                 obj.layerProperties.lengthOfLayers(1, jj)*10^-6 );
%                     ts{jj+1,jj}(:,:,ii) = obj.layerProperties.Ts{ obj.layerProperties.interfaceNum(jj), 2 }(:,:,ii) * ...
%                         expm( -1i*obj.modes.betas{ obj.layerProperties.crossSectionNum(jj+1, 1) }(:,:,ii) * ...
%                                 obj.layerProperties.lengthOfLayers(jj+1)*10^-6 );
%                     % reflect
%                     rs{jj,jj+1}(:,:,ii) = obj.layerProperties.Rs{ obj.layerProperties.interfaceNum(jj),1}(:,:,ii) * ...
%                         expm( -1i*obj.modes.betas{ obj.layerProperties.crossSectionNum(jj,1)}(:,:,ii) * ...
%                                 obj.layerProperties.lengthOfLayers(1,jj)*10^-6 );
%                     rs{jj+1,jj}(:,:,ii) = obj.layerProperties.Rs{ obj.layerProperties.interfaceNum(jj),2}(:,:,ii) * ...
%                         expm( -1i*obj.modes.betas{ obj.layerProperties.crossSectionNum(jj+1,1)}(:,:,ii) * ...
%                                 obj.layerProperties.lengthOfLayers(1,jj+1)*10^-6 );     

                    % DEBUG
%                     fprintf('layer2layer\n');
%                     tic;    
%                     NEWER VERSION, more efficient
%                     transfer
                    ts{jj,jj+1}(:,:,ii) = obj.layerProperties.Ts{ obj.layerProperties.interfaceNum(jj), 1 }(:,:,ii) * ...      
                        diag( exp( -1i*diag( obj.modes.betas{ obj.layerProperties.crossSectionNum(jj) }(:,:,ii) ) .* ...
                                obj.layerProperties.lengthOfLayers(jj)*10^-6 ) );
                    ts{jj+1,jj}(:,:,ii) = obj.layerProperties.Ts{ obj.layerProperties.interfaceNum(jj), 2 }(:,:,ii) * ...
                        diag( exp( -1i*diag( obj.modes.betas{ obj.layerProperties.crossSectionNum(jj+1) }(:,:,ii) ) .* ...
                                obj.layerProperties.lengthOfLayers(jj+1)*10^-6 ) );      
%                     reflect
                    rs{jj,jj+1}(:,:,ii) = obj.layerProperties.Rs{ obj.layerProperties.interfaceNum(jj),1}(:,:,ii) * ...
                        diag( exp( -1i*diag( obj.modes.betas{ obj.layerProperties.crossSectionNum(jj) }(:,:,ii) ) .* ...
                                obj.layerProperties.lengthOfLayers(jj)*10^-6 ) );
                    rs{jj+1,jj}(:,:,ii) = obj.layerProperties.Rs{ obj.layerProperties.interfaceNum(jj),2}(:,:,ii) * ...
                        diag( exp( -1i*diag( obj.modes.betas{ obj.layerProperties.crossSectionNum(jj+1) }(:,:,ii) ) .* ...
                                obj.layerProperties.lengthOfLayers(jj+1)*10^-6 ) );
%                     toc;    % DEBUG
                            
%                     fprintf('inverts\n');
%                     tic;    % DEBUG
                    % recursive formula for calculating inter-layer TS
                    % and RS
                    if jj>1                                                      %Also calculate from front to arbitrary layer for building fields
                        topvar = ts{jj,jj+1}(:,:,ii)/( eye(numodes, numodes) - rs{jj,1}(:,:,ii) * rs{jj,jj+1}(:,:,ii) );
                        ts{1,jj+1}(:,:,ii) = topvar * ts{1,jj}(:,:,ii);
                        rs{jj+1,1}(:,:,ii) = topvar * rs{jj,1}(:,:,ii) * ts{jj+1,jj}(:,:,ii) + rs{jj+1,jj}(:,:,ii);

                        botvar = ts{jj,1}(:,:,ii)/( eye(numodes,numodes) - rs{jj,jj+1}(:,:,ii) * rs{jj,1}(:,:,ii) );
                        rs{1,jj+1}(:,:,ii) = botvar * rs{jj,jj+1}(:,:,ii)*ts{1,jj}(:,:,ii) + rs{1,jj}(:,:,ii);
                        ts{jj+1,1}(:,:,ii) = botvar * ts{jj+1,jj}(:,:,ii);
                    end
%                     toc;    % DEBUG

                end
            end
            obj.scatterProperties.ts = ts;
            obj.scatterProperties.rs = rs;

            toc;    % DEBUG
            fprintf('...done\n\n');
            
            %Get each forward (As) and backward (Bs) coefficients
            As = cell(obj.layerProperties.numOfLayers,1);                       %Instantiate cells
            Bs = cell(obj.layerProperties.numOfLayers,1);
            As{1} = zeros(numodes,length(obj.domain.k0));                       %Initializes coefficients for input field to all zeros (we will add a value to the mode we wish to launch as a source)
            Bs{obj.layerProperties.numOfLayers} = ...
                zeros(numodes,length(obj.domain.k0));                           %No light coming backwards from the end of the grating coming in the back
            norm = zeros(1,length(obj.domain.k0));

            for ii = 1:length(obj.domain.k0)                                    %For each wavelength
                As{1}(obj.sources.modeNumber,ii) = 1; 
                As{obj.layerProperties.numOfLayers}(:,ii) =...
                    ts{1,obj.layerProperties.numOfLayers}(:,:,ii)*As{1}(:,ii);  %Coefficient forward at last layer 
                Bs{1}(:,ii) = rs{1,obj.layerProperties.numOfLayers}...
                    (:,:,ii)*As{1}(:,ii);                                       %Coefficient coming backwards out of the first layer

                %Get coefficients in each layer, Bs (backward) and As (forward)
               for p = obj.layerProperties.numOfInterfaces:-1:2                 %For coefficients at each interface except first
                    As{p}(:,ii) = (eye(numodes,numodes)-rs{p,1}(:,:,ii)*...
                        rs{p,p+1}(:,:,ii))\(ts{1,p}(:,:,ii)*As{1}(:,ii)+...
                        rs{p,1}(:,:,ii)*ts{p+1,p}(:,:,ii)*Bs{p+1}(:,ii));
                    Bs{p}(:,ii) = rs{p,p+1}(:,:,ii)*As{p}(:,ii) + ...
                        ts{p+1,p}(:,:,ii)*Bs{p+1}(:,ii);
                end

            end
            obj.scatterProperties.coeff.As = As;
            obj.scatterProperties.coeff.Bs = Bs;

            % update z vector
            z_vec = obj.domain.z;
%             n_dz    = size(obj.diel, 2);                                % number of disc. in z
%             len_z   = sum( obj.layerProperties.lengthOfLayers(:) );     % total length of z domain
%             z_vec   = linspace( 0, len_z, n_dz );                       % z coordinate vector
%             obj.domain.z = z_vec;
            
            % update the discretization
%             obj.domain.discretization(2) = z_vec(2) - z_vec(1);                              % disc. size in z
            
            %Build fields
%             Fy = zeros(numodes,round(sum(obj.layerProperties.lengthOfLayers)... %Instantiate Fy (which is Ey for TE and Hy for TM)
%                 /obj.domain.discretization(2)),length(obj.domain.k0));
            Fy      = zeros(size(obj.diel,1), size(obj.diel, 2), length(obj.domain.k0) );         % Instantiate field Fy (ey for TE Hy for TM)
            
            for ii = 1:length(obj.domain.k0)                                
                % for each wavelength
%                 zindex = 1; % UNUSED
                
                % starting z position
                z_pos_start = 0;
                
                for jj = 1:obj.layerProperties.numOfLayers
                    % for each layer
                   
                    % grab end position of this layer
                    z_pos_end = z_pos_start + obj.layerProperties.lengthOfLayers(jj);
                    
                    % get current z vector chunk
                    % OLD
%                     z = ([1:round(obj.layerProperties.lengthOfLayers(1,jj)/...
%                         obj.domain.discretization(2))]-0.5)*...
%                         obj.domain.discretization(2);                           %Z-vector in microns
                    % NEW
                    z_indx      = z_vec >= z_pos_start & z_vec <= z_pos_end;           % indices of current z vector
                    diff_z_pos  = z_vec( z_indx ) - z_pos_start;                       % differential z, starting from prev. layer interface
%                     if jj == 1
%                         z_indx      = z_vec >= z_pos_start & z_vec <= z_pos_end;           % indices of current z vector
%                         diff_z_pos  = z_vec( z_indx );                                                       % differential z, starting from prev. layer interface
%                     else
%                         z_indx      = z_vec >= z_pos_start & z_vec <= z_pos_end;
%                         diff_z_pos  = z_vec( z_indx ) - obj.layerProperties.lengthOfLayers(jj);
%                     end
   
%                     % OLD
%                     prop = exp(-1i*10^-6*obj.modes.betas{...
%                         obj.layerProperties.crossSectionNum(jj,1)}(:,:,ii)*...
%                         (diag(z)*ones(length(z),numodes)).');                   %Matrix where row i is e^(-j beta_i z) 
                    % NEW
                    prop = exp(-1i*10^-6*obj.modes.betas{...
                        obj.layerProperties.crossSectionNum(jj,1)}(:,:,ii)*...
                        (repmat(diff_z_pos, [numodes, 1] )) ); % (diag(z)*ones(length(z),numodes)).');                   %Matrix where row i is e^(-j beta_i z) 
                    
                    forwardfields = obj.modes.fields{...                        %Forward fields are propagated forward
                        obj.layerProperties.crossSectionNum(jj,1)}(:,:,ii)*...
                        diag(obj.scatterProperties.coeff.As{jj}(:,ii))*prop;
                    
                    backwardfields = fliplr(obj.modes.fields{...                %Backward fields are propagated forward and then flipped left right for numerical stability
                        obj.layerProperties.crossSectionNum(jj,1)}(:,:,ii)*...
                        diag(obj.scatterProperties.coeff.Bs{jj}(:,ii))*prop);
                    
                    % OLD
%                     Fy(:,zindex:round(sum(...                                   %Places the fields into Fy
%                         obj.layerProperties.lengthOfLayers(1,1:jj))/...
%                         obj.domain.discretization(2)),ii) = forwardfields+...
%                         backwardfields; 
                    % NEW
                    Fy(:, z_indx, ii) = forwardfields + backwardfields;
                    
%                     % OLD
%                     zindex = round(sum(obj.layerProperties.lengthOfLayers...    %Keep track of the zindex where we are in Fy
%                         (1,1:jj))/obj.domain.discretization(2)+1); 
                    
                    % NEW
                    z_pos_start = z_pos_end;    % move to next layer
               end
                
               %Normalize fields to input power here
               
            end
            
            % I'm assuming this is supposed to be a function of # of modes,
            % not length of field
%             % old
%             obj.scatterProperties.PowerIn = zeros(length(obj.domain.k0),length(Fy(:,1)));
%             obj.scatterProperties.PowerRefl = zeros(length(obj.domain.k0),length(Fy(:,1)));
%             obj.scatterProperties.PowerTrans = zeros(length(obj.domain.k0),length(Fy(:,1)));
            % new
            obj.scatterProperties.PowerIn       = zeros( length(obj.domain.k0), numodes );
            obj.scatterProperties.PowerRefl     = zeros( length(obj.domain.k0), numodes );
            obj.scatterProperties.PowerTrans    = zeros( length(obj.domain.k0), numodes );
            
            if obj.domain.polarization == 0                                     %If TE mode 
               obj.fullFields.Ey = Fy;                                          %Store Ey fields
               for ii = 1:length(obj.domain.k0)
                   obj.scatterProperties.PowerIn(ii,:) = ...                        %Power launched in the left side of the domain (row of the power in each mode)
                       0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...               %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                       real(sum(...
                       obj.modes.B{1,1}(:,:,ii)*...
                       abs(obj.modes.fields{1,1}(:,:,ii)).^2*...
                       conj(obj.modes.betas{1,1}(:,:,ii))*...
                       diag(abs(obj.scatterProperties.coeff.As{1,1}...
                       (:,ii)).^2)))*obj.domain.discretization(1)*10^-6;

                   obj.scatterProperties.PowerRefl(ii,:) = ...                      %Power reflected back out the left side of the domain (row of the power in each mode)
                      0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...               %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                      real(sum(...
                      obj.modes.B{1,1}(:,:,ii)*...
                      abs(obj.modes.fields{1,1}(:,:,ii)).^2*...
                      conj(obj.modes.betas{1,1}(:,:,ii))*...
                      diag(abs(diag(diag(exp(-1j.*obj.modes.betas{1,1}...
                      (:,:,1)*obj.layerProperties.lengthOfLayers(1,1)*...
                      10^-6)))*obj.scatterProperties.coeff.Bs{1,1}(:,ii)).^2)...
                       ))*obj.domain.discretization(1)*10^-6; 

                   obj.scatterProperties.PowerTrans(ii,:) = ...                     %Power coming out the right side of the domain (row of the power in each mode)
                      0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...               %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                      real(sum(...
                      obj.modes.B{obj.layerProperties.crossSectionNum(end,1)...
                      ,1}(:,:,ii)*abs(obj.modes.fields{...
                      obj.layerProperties.crossSectionNum(end,1),1}(:,:,ii))...
                      .^2*conj(obj.modes.betas{...
                      obj.layerProperties.crossSectionNum(...
                      end,1),1}(:,:,ii))*diag(abs(diag(diag(exp(-1j.*...
                      obj.modes.betas{obj.layerProperties.crossSectionNum(...
                      end,1),1}(:,:,1)*obj.layerProperties.lengthOfLayers(1,...
                      end)*10^-6)))*obj.scatterProperties.coeff.As{end,1}...
                      (:,ii)).^2)))*obj.domain.discretization(1)*10^-6;

                  %Get powers up and down
                    topLoc = round((obj.domain.domain(1)-...                    %index 2pts below the pmls 
                        (obj.domain.pml(1)+2*obj.domain.discretization(1)))./...
                        obj.domain.discretization(1));
                    botLoc = round((...                                         %index 2pts above the pmls 
                        (obj.domain.pml(1)+2*obj.domain.discretization(1)))./...
                        obj.domain.discretization(1));

                    EyPlaneTop = obj.fullFields.Ey(topLoc,:,ii);                
                    HzPlaneTop = ...
                        1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                        obj.constants.mu0)*...
                        (obj.fullFields.Ey(topLoc+1,:,ii)-...
                        obj.fullFields.Ey(topLoc-1,:,ii))./...
                        (2.*obj.domain.discretization(1)*10^-6);

                    obj.scatterProperties.EyPlaneTop(ii,:) = EyPlaneTop;  
                    obj.scatterProperties.HzPlaneTop(ii,:) = HzPlaneTop;

                    EyPlaneBot = obj.fullFields.Ey(botLoc,:,ii);
                    HzPlaneBot = ...
                        1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                        obj.constants.mu0)*...
                        (obj.fullFields.Ey(botLoc+1,:,ii)-...
                        obj.fullFields.Ey(botLoc-1,:,ii))./...
                        (2.*obj.domain.discretization(1)*10^-6);

                    obj.scatterProperties.EyPlaneBot(ii,:) = EyPlaneBot;
                    obj.scatterProperties.HzPlaneBot(ii,:) = HzPlaneBot;

                    obj.scatterProperties.PowerTop(ii) = ...
                        abs(.5*real(EyPlaneTop*HzPlaneTop'*...
                        obj.domain.discretization(2)*10^-6));

                    obj.scatterProperties.PowerBot(ii) = ...
                        abs(.5*real(EyPlaneBot*HzPlaneBot'*...
                        obj.domain.discretization(2)*10^-6));
               end
                   
            elseif obj.domain.polarization == 1                                 %If TM mode
               obj.fullFields.Hy = Fy;
               
               for ii = 1:length(obj.domain.k0)
                   
                   obj.scatterProperties.PowerIn(ii,:) = ...                    %Power launched in the left side of the domain (row of the power in each mode)
                       0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...          %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                       real(sum(...
                       obj.modes.B{1,1}(:,:,ii)*...
                       abs(obj.modes.fields{1,1}(:,:,ii)).^2*...
                       (obj.modes.betas{1,1}(:,:,ii))*...
                       diag(abs(obj.scatterProperties.coeff.As{1,1}...
                       (:,ii)).^2)))*obj.domain.discretization(1)*10^-6;

                   obj.scatterProperties.PowerRefl(ii,:) = ...                  %Power reflected back out the left side of the domain (row of the power in each mode)
                      0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...           %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                      real(sum(...
                      obj.modes.B{1,1}(:,:,ii)*...
                      abs(obj.modes.fields{1,1}(:,:,ii)).^2*...
                      (obj.modes.betas{1,1}(:,:,ii))*...
                      diag(abs(diag(diag(exp(-1j.*obj.modes.betas{1,1}(:,:,1....
                      )*obj.layerProperties.lengthOfLayers(1,1)*10^-6)))*...
                      obj.scatterProperties.coeff.Bs{1,1}(:,ii)).^2)...
                       ))*obj.domain.discretization(1)*10^-6; 

                   obj.scatterProperties.PowerTrans(ii,:) = ...                 %Power coming out the right side of the domain (row of the power in each mode)
                      0.5./(obj.domain.k0(ii)*10^6*obj.constants.c)*...           %omega in SI = obj.domain.k0(ii)*10^6*obj.constants.c
                      real(sum(...
                      obj.modes.B{obj.layerProperties.crossSectionNum(end,1),...
                      1}(:,:,ii)*abs(obj.modes.fields{...
                      obj.layerProperties.crossSectionNum(end,1),1}(:,:,ii)...
                      ).^2*(obj.modes.betas{...
                      obj.layerProperties.crossSectionNum(...
                      end,1),1}(:,:,ii))*diag(abs(diag(diag(exp(-1j.*...
                      obj.modes.betas{obj.layerProperties.crossSectionNum(end...
                      ,1),1}(:,:,1)*obj.layerProperties.lengthOfLayers(1,end...
                      )*10^-6)))*obj.scatterProperties.coeff.As{end,1}(:,ii)...
                      ).^2)))*obj.domain.discretization(1)*10^-6;
                  
                  %Get powers up and down
                    topLoc = round((obj.domain.domain(1)-...
                        (obj.domain.pml(1)+2*obj.domain.discretization(1)))./...
                        obj.domain.discretization(1));
                    botLoc = round((...
                        (obj.domain.pml(1)+2*obj.domain.discretization(1)))./...
                        obj.domain.discretization(1));

                    HyPlaneTop = obj.fullFields.Hy(topLoc,:,ii);
                    EzPlaneTop = ...
                        -1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                        obj.constants.eps0)*1./(obj.diel(topLoc,:)).^2.*...
                        (obj.fullFields.Hy(topLoc+1,:,ii)-...
                        obj.fullFields.Hy(topLoc-1,:,ii))./...
                        (2.*obj.domain.discretization(1)*10^-6);

                    obj.scatterProperties.HyPlaneTop(ii,:) = HyPlaneTop;
                    obj.scatterProperties.EzPlaneTop(ii,:) = EzPlaneTop;

                    HyPlaneBot = obj.fullFields.Hy(botLoc,:,ii);
                    EzPlaneBot = ...
                        -1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                        obj.constants.eps0)*1./(obj.diel(topLoc,:)).^2.*...
                        (obj.fullFields.Hy(botLoc+1,:,ii)-...
                        obj.fullFields.Hy(botLoc-1,:,ii))./...
                        (2.*obj.domain.discretization(1)*10^-6);

                    obj.scatterProperties.HyPlaneBot(ii,:) = HyPlaneBot;
                    obj.scatterProperties.EzPlaneBot(ii,:) = EzPlaneBot;

                    obj.scatterProperties.PowerTop(ii) = ...
                        abs(.5*real(EzPlaneTop*HyPlaneTop'*...
                        obj.domain.discretization(2)*10^-6));

                    obj.scatterProperties.PowerBot(ii) = ...
                        abs(.5*real(EzPlaneBot*HyPlaneBot'*...
                        obj.domain.discretization(2)*10^-6));
                    
               end
                 
            else
                Error('\n No specified Polarization');                          %This shouldn't happen as it should default to TE
            end
            
            temp1 = size(obj.diel);
            temp2 = size(Fy);
            
            if temp1(2) ~= temp2(2)
                obj = obj.setDiel;
                fprintf(['\nWarning: obj.diel is different size than full'... 
                    ' fields.  Resetting obj.diel']);
            end
              
        end                               %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%    Get field and modes at observation planes     %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = addObsPlane(obj,varargin)
            
            %Works for overlap with only one mode
            %Works only in TE mode so far
            %Works only for vertical observation planes
            %Also doesn't work for curved waveguide modes (not sure if this
            %   will ever be implemented)
            
            fprintf('\nAdding observation plane.\n');
            
            
            obsPlaneFields = {'extents','modesToOverlap',...
                'plotModeToOverlap'};
            
            numObsPlane = length(obj.obsPlanes) + 1;

            for k = 1:2:size(varargin,2)                                        % 
                varName = varargin{k};
                switch(varName)
                    case obsPlaneFields
                        eval(['obj.obsPlanes(' num2str(numObsPlane) ').' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid'...
                            ' parameter "' varargin{k} '".  Check EME'...
                            ' instructions by typing ''help emeSim''.']);
                end
            end
            
            %Ey(x,z)
            extents = obj.obsPlanes(numObsPlane).extents; 
            dx = obj.domain.discretization(1); 
            dz = obj.domain.discretization(2); 
            xLoc1 = round(extents(3)/dx);
            xLoc2 = round(extents(4)/dx);
            zLoc1 = round(extents(1)/dz);
            zLoc2 = round(extents(2)/dz); 
            
            [xSize, zSize] = size(obj.fullFields.Ey(:,:,1));
            
            if xLoc1 > xSize-1
                xLoc1 = xSize-1;
                fprintf(['Warning:  You selected an x-component of your observation plane'...
                    ' in the PMLs or outside the domain.  Setting near edge \n']);
            end
            if xLoc2 > xSize-1
                xLoc2 = xSize-1;
                 fprintf(['Warning:  You selected an x-component of your observation plane'...
                    ' in the PMLs or outside the domain.  Setting near edge \n']);
            end
            if zLoc1 > zSize-1
                zLoc1 = zSize-1;
                fprintf(['Warning:  You selected an z-component of your observation plane'...
                    ' in the PMLs or outside the domain.  Setting near edge \n']);
            end
            if zLoc2 >zSize-1
                zLoc2 = zSize-1;
                fprintf(['Warning:  You selected an z-component of your observation plane'...
                    ' in the PMLs or outside the domain.  Setting near edge \n']);
            end
            
            
            for ii = 1:length(obj.domain.k0)
                if xLoc1 == xLoc2
                    obj.obsPlanes(numObsPlane).field(:,ii) = obj.fullFields.Ey(xLoc1,zLoc1:zLoc2,ii).'; 
                    n = obj.diel(xLoc1,zLoc1:zLoc2);
                    %n(1:zLoc1) = n(zLoc1); n(zLoc2:end) = n(zLoc2);
                    obj.obsPlanes(numObsPlane).n = n;
                    dd = dz;
                    axisString = 'z';

                elseif zLoc1 == zLoc2 
                    obj.obsPlanes(numObsPlane).field(:,ii) = obj.fullFields.Ey(xLoc1:xLoc2,zLoc1,ii);
                    n = obj.diel(xLoc1:xLoc2,zLoc1);
                    %n(1:xLoc1) = n(xLoc1); n(xLoc2:end) = n(xLoc2);
                    obj.obsPlanes(numObsPlane).n = n;
                    dd = dx;
                    axisString = 'x';
                else
                    error('\nObservation plane extents error \n');
                end


                theMode = modeslv(n, dd, obj.domain.k0(ii), obj.domain.polarization);
                obj.obsPlanes(numObsPlane).mode(:,ii) = theMode(:,obj.obsPlanes(numObsPlane).modesToOverlap);


                %Get total power in that cross-section
                        %Get total power in that cross-section

                Efield = (obj.obsPlanes(numObsPlane).field(:,ii)).';
                if strcmp(axisString, 'z')       %This is if the obsPlane is horizontal
                    error('\n obsPlanes does not support horizontal obs planes yet\n');
                    Hfield =  (1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                            obj.constants.mu0)*...
                            (obj.fullFields.Ey(xLoc1+1,zLoc1:zLoc2,ii)-...
                            obj.fullFields.Ey(xLoc1-1,zLoc1:zLoc2,ii))./...
                            (2.*obj.domain.discretization(1)*10^-6)).';
                elseif strcmp(axisString, 'x')     %This is if the obsPlane is vertical
                    Hfield =  (1j.*1/(obj.domain.k0(ii)*10^6*obj.constants.c*...
                            obj.constants.mu0)*...
                            (obj.fullFields.Ey(xLoc1:xLoc2,zLoc1+1,ii)-...
                            obj.fullFields.Ey(xLoc1:xLoc2,zLoc1-1,ii))./...
                            (2.*obj.domain.discretization(1)*10^-6)).';
                else
                    error('\nMode overlap does not support angled obs planes yet\n');
                end

                Power(ii) = ...
                        abs(.5*real(Efield*Hfield'*...
                        obj.domain.discretization(2)*10^-6));


                Emode = (obj.obsPlanes(numObsPlane).mode(:,ii)).';    
                obj.obsPlanes(numObsPlane).powInMode(ii) = ...
                     Power(ii)*abs(Efield*Emode')^2/...
                                    abs(real((Efield*Efield')*(Emode*Emode')));


                fprintf(['\nPower in obs plane is ' num2str(Power(ii)) '\n']);    
            
            end
            if strcmp(obj.obsPlanes(numObsPlane).plotModeToOverlap, 'yes')
                xx = (xLoc1:1:xLoc2)*dd;
                figure; plot(xx,abs(obj.obsPlanes(numObsPlane).mode),'LineWidth',2);
                hold on;
                plot(xx,abs(obj.obsPlanes(numObsPlane).field),'LineWidth',2);
                hold off;
                xlabel([axisString ' (\mum)']);
                ylabel('|Field| (a.u.)'); 
                legend('Mode','Field');
                title(['Field and Mode to Overlap at Obs Plane ' num2str(numObsPlane)]);
                set(gca,'FontSize', 14);
                set(gcf,'Color','w'); grid on;
            end
            
        end                             %-------%%%
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%             Fiber Mode Overlap                   %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = fiberOverlap(obj,varargin)
            
            fprintf('\nFiber overlap.\n');

            overlapFields = {'overlapDir','MFD','angleVec','zOffset','nClad'};          %
            for k = 1:2:size(varargin,2)                                        % 
                varName = varargin{k};
                switch(varName)
                    case overlapFields
                        eval(['obj.fiberCoup.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid'...
                            ' parameter "' varargin{k} '".  Check EME'...
                            ' instructions by typing ''help emeSim''.']);
                end
            end
            if ~isfield(obj.fiberCoup,'overlapDir')                             %If modeNumber isn't set, set to 1, which is usually the fundamental mode
                obj.fiberCoup.overlapDir = ...
                    'up';                                                       %Default to up          
                fprintf(['overlapDir was not set.  Setting to default of up'...
                    ' \n']);
            end
            if ~isfield(obj.fiberCoup,'angleVec')                               %If angleVec isn't set
                obj.fiberCoup.angleVec = -30:.5:30;                             %Default to            
                fprintf(['angleVec was not set.  Setting to default of'...
                    ' -30 to 30 \n']);
            end
            if ~isfield(obj.fiberCoup,'zOffset')                                %If zOffset isn't set 
                obj.fiberCoup.zOffset = -5:.5:25;                               %Default to          
                fprintf(['zOffset was not set.  Setting to default of'...
                    ' -5 to 25 \n']);
            end
            if ~isfield(obj.fiberCoup,'nClad')                                  %If nClad isn't set 
                obj.fiberCoup.nClad = 1;                                        %Default to air       
                fprintf(['nClad was not set.  Setting to default of'...
                    ' 1 (i.e., air) \n']);
            end
            
            if strcmp(obj.fiberCoup.overlapDir,'up')
                fprintf('\nOverlapping for light going up.\n');
            elseif strcmp(obj.fiberCoup.overlapDir, 'down')
                fprintf('\nOverlapping for light going down.\n');
            else
                error('Invalid input for overlap. Must be "up" or "down"');
            end
            if ~isfield(obj.fiberCoup,'MFD')                                    %If MFD isn't set, set it to 10.4.
                obj.fiberCoup.MFD = 10.4;
                fprintf('MFD was not set. Setting to 10.4 microns.\n');
            end
            
            angleVec = obj.fiberCoup.angleVec;
            zOff = obj.fiberCoup.zOffset;
            
            % OLD
%             z = obj.domain.discretization(2):obj.domain.discretization(2):...
%                 obj.domain.domain(2);
            % NEW
            z = obj.domain.z;
            
            for  ii = 1:length(obj.domain.k0)
                if obj.domain.polarization == 0
                    if strcmp(obj.fiberCoup.overlapDir, 'up')  
                        Ey = obj.scatterProperties.EyPlaneTop(ii,:);
                        Hz = obj.scatterProperties.HzPlaneTop(ii,:);  
                        P2 = obj.scatterProperties.PowerTop(ii);
                        plusOrMinus = -1;
                    elseif strcmp(obj.fiberCoup.overlapDir, 'down') 
                        Ey = obj.scatterProperties.EyPlaneBot(ii,:);
                        Hz = obj.scatterProperties.HzPlaneBot(ii,:);
                        P2 = obj.scatterProperties.PowerBot(ii);
                        plusOrMinus = 1;
                    end
                    for jj = 1:length(zOff)
                        for kk = 1:length(angleVec)
                            F = fiberModeGaussian(obj.fiberCoup.MFD/2.,...
                                obj.domain.k0(ii)*obj.fiberCoup.nClad...
                                ,0,0,z-zOff(jj),angleVec(kk),0);
                            ey = F.Ey; hz = F.Hz;
                            
                            %{
                            kxfib = obj.domain.k0(ii)*obj.fiberCoup.nClad * sin(angleVec(kk)/180*pi);
                            temp = exp(-1j*z*kxfib);   
                            [ey, hz, ~] = gaussbeam1d(obj.fiberCoup.MFD/2, obj.domain.k0(ii), obj.fiberCoup.nClad, z-zOff(jj));
                            ey = ey.*temp;
                            hz = hz.*temp;
                            plusOrMinus = -plusOrMinus;
                            %}
                            
                                                                                %Scalar version: obj.fiberCoup.coupIn2(jj,kk,ii) = P2*abs(Ey*ey')^2/((ey*ey')*(Ey*Ey'));
                            obj.fiberCoup.coup(jj,kk,ii) = ...
                                P2*abs(Ey*hz'+plusOrMinus*Hz*ey')^2/...
                                abs(real((ey*hz'+hz*ey')*(Ey*Hz'+Hz*Ey')));
                        end
                    end

                elseif obj.domain.polarization == 1
                    if strcmp(obj.fiberCoup.overlapDir, 'up')
                        Hy = obj.scatterProperties.HyPlaneTop(ii,:);
                        Ez = obj.scatterProperties.EzPlaneTop(ii,:);  
                        P2 = obj.scatterProperties.PowerTop(ii);
                        plusOrMinus = -1;
                    elseif strcmp(obj.fiberCoup.overlapDir, 'down')
                        Hy = obj.scatterProperties.HyPlaneBot(ii,:);
                        Ez = obj.scatterProperties.EzPlaneBot(ii,:);
                        P2 = obj.scatterProperties.PowerBot(ii);
                        plusOrMinus = 1;
                    end
                    for jj = 1:length(zOff)
                        for kk = 1:length(angleVec)
                            F = fiberModeGaussian(obj.fiberCoup.MFD/2.,...
                                obj.domain.k0(ii)*obj.fiberCoup.nClad...
                                ,1,0,z-zOff(jj),angleVec(kk),0);
                            hy = F.Hy; ez = F.Ez; 
                                                                                %Scalar version: obj.fiberCoup.coupIn2(jj,kk,ii) = P2*abs(Hy*hy')^2/((hy*hy')*(Hy*Hy')); 
                            obj.fiberCoup.coup(jj,kk,ii) = ...
                                P2*abs(Hy*ez'+plusOrMinus*Ez*hy')^2/...
                                abs(real((hy*ez'+ez*hy')*(Hy*Ez'+Ez*Hy')));
                        end
                    end
                end
                
            end
            
            centerLambdaIndex = round(length(obj.domain.k0)/2);
            % OLD CODE
%             [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); 
%             [obj.fiberCoup.optCoup, I2] = max(tempmax);
%             obj.fiberCoup.optZOffsetIndex = I1(I2);
%             obj.fiberCoup.optAngleIndex = I2;
%             obj.fiberCoup.optZOffset = obj.fiberCoup.zOffset(I1(I2));
%             obj.fiberCoup.optAngle = obj.fiberCoup.angleVec(I2);
            
            % NEW CODE
            coupling_center_lambda = obj.fiberCoup.coup(:,:,centerLambdaIndex); % dimensions zoffset vs. angle
            [ obj.fiberCoup.optCoup, indx_opt ] = max(coupling_center_lambda(:)); 
            [ indx_zoff, indx_angle ] = ind2sub( size(coupling_center_lambda), indx_opt );
            obj.fiberCoup.optZOffset = obj.fiberCoup.zOffset( indx_zoff );
            obj.fiberCoup.optAngle = obj.fiberCoup.angleVec( indx_angle );
 
        end                            %-------%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%    Run Multiple function for Full Simulation     %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = runSimulation(obj,varargin)
            % Runs the full scatter/propagate simulation
            %
            % Inputs
            %   varargin
            %       Desc: name value pairs. All are optional
            %       'plotSource'
            %           'yes' or 'no' to plot the source mode
            %       'n_modes'
            %           number of modes to solve for. Defaults to all modes
            
            runSimFields    = {'plotSource', 'n_modes'};   % valid optional fieldnames
            OPTS            = struct;
            for k = 1:2:size(varargin,2)                                        
                varName = varargin{k};
                switch(varName)
                    case runSimFields
%                         eval([varName ' = varargin{k+1};']);
                        OPTS = setfield( OPTS, varName, varargin{k+1} );
                    otherwise
                        error(['Trying to assign value to invalid'...
                            ' parameter "' varargin{k} '".  Check EME'...
                            ' instructions by typing ''help emeSim''.']);
                end
            end
            % set defaults
            if ~isfield( OPTS, 'plotSource' )
                OPTS.plotSource = 'yes';
            end
            
            % Convert dielectric
            tic;    % DEBUG
            fprintf('Converting dielectric...\n');
            obj = obj.convertDiel();   
            toc;    % DEBUG
            fprintf('Done\n\n');

            
            % Calculate mode/mode overlaps
            tic;    % DEBUG
            fprintf('Calculating modes...\n');
            if ~isfield( OPTS, 'n_modes' )
                obj = obj.getModes(); 
            else
                obj = obj.getModes( OPTS.n_modes ); 
            end
            toc;    % DEBUG
            fprintf('Done\n\n');
            
            
            % Add source
            tic;    % DEBUG
            fprintf('Adding source...\n');
            obj = obj.addSource('modeNumber', 1, 'plotSource', OPTS.plotSource);  
            toc;    % DEBUG
            fprintf('Done\n\n');
            
            % Run scatter
            tic;    % DEBUG
            fprintf('Performing scatter\n');
            obj = obj.scatter();    
            toc;    % DEBUG
            fprintf('Done\n\n');
            
        end                         %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%      Set the dielectric via layer properties     %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = setDiel(obj)
       
            dx = obj.domain.discretization(1);
            dz = obj.domain.discretization(2);
            zlength = sum(obj.layerProperties.lengthOfLayers);
            zSize = round(1e8*zlength/dz)/(1e8);
            xSize = length(obj.layerProperties.uniqueCrossSections(:,1));
           
            for ii = 1:obj.layerProperties.numOfLayers
                tempN = obj.layerProperties.uniqueCrossSections(:, ...
                    obj.layerProperties.crossSectionNum(ii));
                
                layerLength = round(1e8*obj.layerProperties.lengthOfLayers(ii)/dz)/(1e8);
                [na Nsect] = meshgrid(1:layerLength,tempN);
                
                if ii == 1
                    n = Nsect;
                else
                    n = [n Nsect];
                end
            end
            
            obj.diel = n;
            obj.domain.domain = [xSize*dx zSize*dz];
              
        end                               %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%     Plots the full field at center wavlength     %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = plotFields(obj)
            %Get the max field value to set caxis
            
            if length(obj.domain.k0)>1
                kLoc = ceil(length(obj.domain.k0)/2.);
            else
                kLoc = 1;
            end
            
            if obj.domain.polarization == 0
                field = real(obj.fullFields.Ey(:,:,kLoc)); 
            elseif pol == 1
                field = real(obj.fullFields.Hy(:,:,kLoc)); 
            else
                error('No polarization');
            end
            
            dz = obj.domain.discretization(2);
            dx = obj.domain.discretization(1);
            
            [xsize, zsize] = size(field);
            
            figure; 
            imagesc(dz*(1:zsize), dx*(1:xsize),field); colormap('redbluedark');
            xlabel('z (\mum )'); ylabel('x (\mum)');
            set(gca,'Ydir','normal'); 
            maxx = max(max(abs(field)));
            caxis([-maxx maxx]);
            title('Re\{Field\}','FontSize', 20); set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
        end                            %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%          Plots the index distribution            %-------%%%
        %%%------------------------------------------------------------------%%%
        function obj = plotDiel(obj)
            %Get the max field value to set caxis
           
            [temp1, temp2] = size(obj.diel);
            
            if temp1 == 0 || temp2 == 0
                obj = obj.setDiel();
            end
            
            dz = obj.domain.discretization(2);
            dx = obj.domain.discretization(1);
            
            [xsize, zsize] = size(obj.diel);
            
            figure; 
            imagesc(dz*(1:zsize), dx*(1:xsize),abs(obj.diel)); 
            xlabel('z (\mum )'); ylabel('x (\mum)');
            set(gca,'Ydir','normal'); 
            title('Index','FontSize', 20); set(gca,'FontSize', 14); set(gcf, 'Color', [1 1 1]); 
            colorbar;
        end                              %-------%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%-------%     Returns vector of coupling vs wavelength     %-------%%%
        %%%------------------------------------------------------------------%%%
        %%% Added by Bohan Zhang
        function [lambdaAll, coupVLambda, optAngle] = couplingVLambda(obj)
            % returns the wavelengths, best coupling efficiency vs wavelength, and best
            % coupling angles
           
            % get wavelengths simulated
            lambda1 = obj.domain.wavelengthSpectrum(1);
            lambda2 = obj.domain.wavelengthSpectrum(2);
            dlambda = obj.domain.wavelengthSpectrum(3);
    
            centerLambdaIndex = round(length(obj.domain.k0)/2);
            
            lambdaAll = lambda1:dlambda:lambda2;
            
            % if multiple wavelengths simulated?
%             if lambda1 ~= lambda2
                [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); 
                [maxCoup, I2] = max(tempmax);
                optZOffset = obj.fiberCoup.zOffset(I1(I2));
                optAngle = obj.fiberCoup.angleVec(I2);

                % Get the efficiency at that z-offset and angle for each wavelength
                for jj = 1:length(lambda1:dlambda:lambda2)
                    coupVLambda(jj)= obj.fiberCoup.coup(I1(I2),I2,jj);
                end
% 
%             plot(lambda1:dlambda:lambda2,coupVLambda, 'LineWidth',1.5); xlabel('Wavelength (\mum)','FontSize', 14);
%             ylabel('Coupling Efficiency','FontSize', 14);  axis([lambda1 lambda2 0 1]);
%             set(gcf,'Color', 'w'); set(gca,'FontSize',12); title('Grating Response','FontSize',14); grid on;
%             else
%                 
%             fprintf('\n You only simulated one wavelength \n');
%             end
            
        end
        
        % -----------------------------------------------------------------
        % Getter function for retrieving layer lengths
        % -----------------------------------------------------------------
        function [ layerLengths ] = getLayerLengths(obj)
            layerLengths = obj.layerProperties.lengthOfLayers;
        end
        
        
        % -----------------------------------------------------------------
        % Setter function for setting layer lengths
        % -----------------------------------------------------------------
        function obj = setLayerLengths(obj, layerLengths)
            % Sets the layer lengths
            %
            % Inputs:
            %   layerLengths
            %       type: float, array
            %       desc: length of each layer in um. The size of the array
            %               must be the same as the current # of layers
            
            % set layer lengths if the number of layers is the same
            if length( obj.layerProperties.lengthOfLayers ) == length( layerLengths )
                obj.layerProperties.lengthOfLayers = layerLengths;
                
                % update z vector
                n_dz    = size(obj.diel, 2);                                % number of disc. in z
                len_z   = sum( obj.layerProperties.lengthOfLayers(:) );     % total length of z domain
                z_vec   = linspace( 0, len_z, n_dz+1 );                       % z coordinate vector
                z_vec   = z_vec(2:end);
                obj.domain.z = z_vec;

                % update the discretization
                obj.domain.discretization(2) = z_vec(2) - z_vec(1);                              % disc. size in z
                
                % update the final z
                obj.domain.domain(2) = obj.domain.z(end);
                
                % re-draw the dielectric profile
                start_layer_z   = obj.domain.z(1);
                end_layer_z     = 0;
                x_section_nums  = obj.layerProperties.crossSectionNum;  % alias
                
                for ii = 1:length(x_section_nums)
                   % for each layer cross-section
                   
                   % end position of this layer
                   end_layer_z = end_layer_z + obj.layerProperties.lengthOfLayers(ii);
                   
                   % assign the dielectric
                   obj.diel( :, z_vec >= start_layer_z & z_vec <= end_layer_z ) = ...
                       repmat( obj.layerProperties.uniqueCrossSections(:, x_section_nums(ii)), 1, sum( z_vec >= start_layer_z & z_vec <= end_layer_z ) );
                   
                   % move to next layer
                   start_layer_z    = end_layer_z;
                   
                   
                end
                
            else
                error('In setLayerLengths(), Layer lengths have changed. Layers will not be set');
            end
            
        end
        
    end%methods
       
end%classdef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- These are extra functions that the class above will call --------$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, beta, B, fundbeta, q, n] = modeslv(n, dx, k0, pol, OPTS)
    % 1-D Modesolver
    % Authors: Cale Gentry, bohan
    %
    %   Inputs:
    %       n:   vector of refractive index at each grid point
    %       dx:  discretization of solver (microns)
    %       k0:  2*pi/(wavelength of light in free space) [rad/micron]            
    %       pol: polarization. 0 for TE 1 for TM (TE is default if pol
    %            undefined)
    %       OPTS:
    %           type: struct
    %           desc: Optional inputs in a struct, with these fields:
    %               PML: 
    %                   PML.d: thickness of PML in micron
    %                       DEFAULT = 0.1 
    %                   PML.m: power strength minus 1 (0:linear, 1:quadratic, ...)
    %                       DEFAULT = 1                           
    %                   PML.xim: maximum distance x goes imaginary (at edges of
    %                       computational domain)
    %                       DEFAULT = 20.
    %               n_modes:
    %                   Integer that chooses # of modes to solve for.
    %                   Default is to solve for all modes
    %     
    %   Outputs:
    %       F: field distribution for each mode (each mode is a column)
    %       beta: matrix of propagation constants for each mode (the diagonal)
    %       B: diagonal matrix of PML and material properties used in overlaps
    %       fundbeta: propagation constant of fundamental mode
    %       q: column/row of beta with smallest imaginary component
    %       n: output of nbuild (probably not used much)

    %OPTS
    if nargin < 5
        OPTS = struct;
    end
    
    if ~isfield(OPTS,'average') 
        OPTS.average = 'yes'; 
    end
    
    if ~strcmp(OPTS.average,'no')
        if ~strcmp(OPTS.average,'yes')
            error('OPTS.average must be yes or no');
        end
    end
    
    %Constants
    c = 299792458;                                                              %SI units (m/s)
    mu0 = 4*pi*10^-7;                                                           %Permeability of free space (SI units) [m kg/(A^2 s^2)]
    eps0 = 1/(c^2*mu0);                                                         %Permittivity of free space (SI units) [s^4 A^2/(m^3 kg)]
   
    lambdaum = 2*pi/k0;
    lambda = lambdaum*10^-6;                                                    %Convert wavelength (microns) to SI units (m)
    f = c/lambda;                                                               %frequency of light (1/s)
    w = 2*pi*f;                                                                 %angular frequency of light (rad/s)
    
    TE = 0;                                                                     %Initialize polarization values
    TM = 1;

    %Default Polarization is TE
    if nargin < 4
        pol = TE; 
    end             

    if ~isfield(OPTS,'PML') OPTS.PML = struct; end                              % Make OPTS.PML structure     
    if ~isfield(OPTS.PML,'d') OPTS.PML.d = .1; end                              % Default PML thickness (micron)
    if ~isfield(OPTS.PML,'m') OPTS.PML.m = 1; end                               % Default power strength (0:line, 1:quadratic...)
    if ~isfield(OPTS.PML,'xim') OPTS.PML.xim = 20.; end                         % Default max imaginary x thickness
    
    dxum = dx;                                                                  % Discretization in microns
    dx = dxum*10^-6;                                                            % Discretization in meters (SI)
    
    N = length(n);                                                              % Number of discretization points in index vector (same length as eps and Fy)
    
    %Material properties
    mu = ones(N-1,1)*mu0;                                                       % Permeability is one cell smaller b/c its associated with (i + 1/2) steps (like H-field)
    eps = n.*n*eps0;                                                            % Permittivity is length N b/c it is associated with i steps (like E-field)
    
    %Averaging
    if strcmp(OPTS.average, 'yes')                                              % Averages adjacent (i + 1/2) steps of permeability to interpolate i step permeabilities
        muavg = ([mu(1);mu] + [mu; mu(end)])*0.5;                               %   The first and last cell of mu vector is averaged with itself
        epsavg = (eps(1:end-1)+eps(2:end))*0.5;
        n = (n(1:end-1)+n(2:end))*.5;
    else
        muavg = ones(N,1)*mu0;                                                  % Else just make mu0 the correct size fromm the start
        epsavg = eps;
    end

    %Create the S variable complex coordinate stretching
    s = ones(N,1);                                                              % Initialize s parameter
    shalf = ones(N-1,1);                                                        % Initialize s parameter on the half grid
    if OPTS.PML.d ~=0                                                           % If 0 distance of PML wasn't requested (ie no PMLs)
        loca = round(OPTS.PML.d/dxum);                                          % Round to the number of discretization points
        s(end-loca:end,1) = 1 - 1i*OPTS.PML.xim/(loca)^OPTS.PML.m*...           % Do bottom section of PML
            ([0:loca]).^OPTS.PML.m;
        s(1:loca+1,1) = (fliplr((s(end-loca:end,1)).')).';                      % Do top section of PML
        shalf(end-loca+1:end,1) = 1 - 1i*OPTS.PML.xim/(loca)^OPTS.PML.m*...     % Half step on bottom section
            ((.5+[0:loca-1])).^OPTS.PML.m;
        shalf(1:loca,1) = (fliplr((shalf(end-loca+1:end,1)).')).';              % Half step on top section
    end

    %Builds tridiagonal matrix (M) to take eigenvalues/vectors of
    if pol == TE
        %Initialize right
        right = muavg(1:end-1)./(s(1:end-1).*shalf.*mu.*dx^2);
        %Initialize left
        left = muavg(2:end)./(s(2:end).*shalf.*mu.*dx^2);
        %Initialize middle
        mid = zeros(N,1);
        %Do all the inner ones
        mid(2:end-1) = muavg(2:end-1).*(w^2.*eps(2:end-1)-1./(s(2:end-1).*dx^2).*(1./(shalf(1:end-1).*mu(1:end-1))+1./(shalf(2:end).*mu(2:end))));
        %Do the corners
        mid(1) = muavg(1).*(w^2.*eps(1)-1/(s(1)*dx^2).* (1./(shalf(1)*mu(1))+1./(shalf(1)*mu(1))));
        mid(end) = muavg(end).*(w^2.*eps(end)-1/(s(end)*dx^2).* (1./(shalf(end)*mu(end))+1./(shalf(end)*mu(end))));
    end
    if pol == TM
        %Initialize right
        right = epsavg(1:end-1)./(s(2:end-1).*shalf(1:end-1).*eps(2:end-1).*dx^2);
        %Initialize left
        left = epsavg(2:end)./(s(2:end-1).*shalf(2:end).*eps(2:end-1).*dx^2);
        %Initialize mid
        mid = epsavg.*(w^2.*mu-1./(shalf.*dx^2).*(1./(s(1:end-1).*eps(1:end-1))+1./(s(2:end).*eps(2:end))));
    end    
    %Build the M matrix
    M = diag(right,1) + diag(left,-1) + diag(mid);

    %Find the eigenvalues/vectors
    %F is the eigenvectors
    %b is the eigenvalues
    if isfield( OPTS, 'n_modes' )
        % Number of modes to use has been specified, run eigs
        [F, b] = eigs(M, OPTS.n_modes, 'sm');
    else
        % Number of modes not specified, default to finding all modes
        [F, b] = eig(M);
    end

    %The square root of the eigenvalues are the propagation constants (beta)
    beta = sqrt(b);
    reals = real(beta);         
    imags = abs(imag(beta));
    beta = reals - 1i*imags;                                                    % This is where I reverse the sign of the imaginary beta

    if pol == TE
        B = diag(s./muavg);
    end
    if pol == TM
        B = diag(shalf./epsavg);
    end
    

    %sort modes from smallest imaginary beta to largest
    betav = diag(beta,0);                                                       % Makes a vertical vector of betas
    imagbeta = abs(imag(betav));                                                % Takes the magnitude of the imaginary 
    matrix = horzcat(imagbeta,betav,F.');                                       % Make the left column the values to be sorted, second column the betas, the rest the fields horizontal
    sorted = sortrows(matrix);                                                  % Sorts the matrix from smallest value of first column
    beta = diag(sorted(:,2));                                                   % Get the betas back in order
    F = (sorted(:,3:end)).';                                                    % Get the fields back 
    q = 1;
    fundbeta = beta(1,1);
    
    %overlap(F1, F2, beta, B)
    %   ((F1).')*B*(F2)*beta;
    F = F*diag(1./sqrt(.5*1/w*diag(F.'*B*F*beta*dx)));                                    %Normalization
    
end






function [n] = nbuild(nlyrs, dlyrsx, dx)
    %1-D Index vector buider
    %   Inputs:
    %       nlyrs: vector of indecies of refraction for each layer
    %       dlyrsx: vector of thicknesses (microns) of each layer
    %       dx:  discretization of solver (microns)


    j = 1;          %Initialize j variable to 1
    last = 0;       %Initialize last variable to 0

    for i = 1:length(nlyrs)     %for each layer 

        last = last + round(dlyrsx(i)/dx);   %make last an integer of previous last plus number of discretizations in current layer
        n(j: last-1) = nlyrs(i);             %make all the values in that layer the same index except the end

        %make the last one in layer averaged with next layer except if it is at
        %the end then there is obviously no layer to average with
        if i~=length(nlyrs)
           n(last) = sqrt((nlyrs(i).^2+nlyrs(i+1).^2)*0.5);   %Average the eps
        end

        j = last +1;
    end

    n = n.'; %Transpose it into a vertical vector
end




