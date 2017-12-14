classdef c_synthesizeGrating
    %C_SYNTHESIZEGRATING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        inputs
        dataObj
        parse
        extract
        synth
        dims
    end
    
    methods
        
        function obj = c_synthesizeGrating(varargin)
            domainFields = {'fileName','optimalAngle'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case domainFields
                        eval(['obj.inputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
        end
        
        function obj = loadData(obj)
            tmp = load(obj.inputs.fileName);
            obj.dataObj = tmp.gratObj;
        end
        
        function obj = parseData(obj)
            % get some variables
            fillVec = obj.dataObj.fillVector;
            ratioVec = obj.dataObj.ratioVector;
            offsetVec = obj.dataObj.offsetVector;
            periodVec = obj.dataObj.periodVector;
            angleMatrix = obj.dataObj.angleMatrix;
            directivityMatrix = obj.dataObj.directivityMatrix;
            scatteringStrengthMatrix = obj.dataObj.scatteringStrengthMatrix;
            optimalAngle = obj.inputs.optimalAngle;
            
            % parse through data and extract matrices
            for II=1:length(fillVec)
                for JJ=1:length(ratioVec)
                    % for each offset (row), find period that gives angle closest to
                    % optimal angle
                    tempAngleMatrix = squeeze(angleMatrix(II,JJ,:,:));
                    [minErrorVec, minAngIdx] = min(abs(tempAngleMatrix-optimalAngle),[],2);
                    dMatrix = squeeze(directivityMatrix(II,JJ,:,:));
                    
                    for KK=1:length(minAngIdx)
                        dVec(KK) = dMatrix(KK,minAngIdx(KK));
                    end
                    
                    [minD, minDIdx] = min(abs(dVec));
                    
                    % keep these variables
                    obj.parse.offsetMatrix(II,JJ) = offsetVec(minDIdx);
                    obj.parse.periodMatrix(II,JJ) = periodVec(minAngIdx(minDIdx));
                    obj.parse.angleMatrix(II,JJ) = tempAngleMatrix(minDIdx,minAngIdx(minDIdx));
                    obj.parse.directivityMatrix(II,JJ) = minD;
                    obj.parse.scatteringStrengthMatrix(II,JJ) = scatteringStrengthMatrix(II,JJ,minDIdx,minAngIdx(minDIdx));                
                end
            end
        end
        
        function obj = extractContours(obj)
            % this method extracts the contours from the simulation plots
            
            %extract some variables
            fillVec = obj.dataObj.fillVector;
            ratioVec = obj.dataObj.ratioVector;
            directivityMatrix = obj.parse.directivityMatrix;
            periodMatrix = obj.parse.periodMatrix; 
            angleMatrix = obj.parse.angleMatrix; 
            scatteringStrengthMatrix = obj.parse.scatteringStrengthMatrix;
            offsetMatrix = obj.parse.offsetMatrix; 
            
            peakLocs = nan(length(fillVec),4);
            logDirectivityDown = log10(abs(1./directivityMatrix));
            myFilter = fspecial('gaussian',[3 3], 1.0);
            filteredLogDirectivityDown = imfilter(logDirectivityDown, myFilter, 'replicate');
            
            for II=1:length(fillVec)
                dVecLine = filteredLogDirectivityDown(II,:);
                if max(dVecLine)>1 %enforce minimum global peak value
                    dVecLine = dVecLine./max(dVecLine);
                else
                    dVecLine = dVecLine*0;
                end
                [~,locs]=findpeaks(dVecLine,ratioVec,'NPeaks',4,'MinPeakHeight',0.5,'MinPeakProminence',0.1);
                if length(locs)>1
                    peakLocs(II,1:length(locs))=locs;
                elseif isempty(locs)
                    peakLocs(II,1:length(locs))=nan;
                else
                    peakLocs(II,1)=locs;
                end
            end
            
            % parse through the data to extract points to fit
            peakLocX_1 = cell(length(ratioVec),size(peakLocs,2)); %cell array
            peakLocY_1 = cell(length(ratioVec),size(peakLocs,2)); %cell array
            
            for II=1:size(peakLocs,2)
                tempVec = peakLocs(:,II);
                [sortedVec sortedIdx] = sort(tempVec);
                tempPeakLocX(:,II) = sortedVec;
                tempPeakLocY(:,II) = fillVec(sortedIdx);
                for JJ=1:size(peakLocX_1,1)
                    idx = find(tempPeakLocX(:,II)==ratioVec(JJ));
                    peakLocX_1{JJ,II}=ratioVec(JJ);
                    if isempty(idx)
                        peakLocY_1{JJ,II}=nan;
                    else
                        peakLocY_1{JJ,II}=tempPeakLocY(idx,II);                        
                    end
                end
            end
            
            % scan over ratioVec (x)
            peakLocs = nan(length(ratioVec),4);
            for II=1:length(ratioVec)
                dVecLine = filteredLogDirectivityDown(:,II);
                if max(dVecLine)>1 %enforce minimum global peak value
                    dVecLine = dVecLine./max(dVecLine);
                else
                    dVecLine = dVecLine*0;
                end
                [~,locs]=findpeaks(dVecLine,fillVec,'NPeaks',4,'MinPeakHeight',0.5,'MinPeakProminence',0.1);
                if length(locs)>1
                    peakLocs(II,1:length(locs))=locs;
                elseif isempty(locs)
                    peakLocs(II,1:length(locs))=nan;
                else
                    peakLocs(II,1)=locs;
                end
            end
            
            % parse through the data to extract points to fit
            peakLocX_2 = cell(length(ratioVec),size(peakLocs,2)); %cell array
            peakLocY_2 = cell(length(ratioVec),size(peakLocs,2)); %cell array
            
            for II=1:size(peakLocs,2)
                for JJ=1:size(peakLocs,1)
                    peakLocX_2{JJ,II}=ratioVec(JJ);
                    peakLocY_2{JJ,II}=peakLocs(JJ,II);
                end
            end
            
            %combine peakLocXY_1 and _2
            for II=1:length(ratioVec)
                combinedPeakLocs_X{II} = ratioVec(II);
                combinedPeakLocs_Y{II} = [];
                for JJ=1:length({peakLocY_1{II,:}})
                    if ~isnan(peakLocY_1{II,JJ})
                        combinedPeakLocs_Y{II} = [combinedPeakLocs_Y{II} peakLocY_1{II,JJ}.'];
                    end
                end
                for JJ=1:length({peakLocY_2{II,:}})
                    if ~isnan(peakLocY_2{II,JJ})
                        combinedPeakLocs_Y{II} = [combinedPeakLocs_Y{II} peakLocY_2{II,JJ}.'];
                    end
                end
            end
            
            % now sort through the simulation results to extract the inverted design
            for II=1:length(combinedPeakLocs_X)
                if II==1
                    yKeep(II)=mean(combinedPeakLocs_Y{II});
                    xKeep(II)=combinedPeakLocs_X{II};
                    yTempKeep=yKeep(II);
                else
                    yTemp = combinedPeakLocs_Y{II};
                    for JJ=1:length(yTemp)
                        
                        if mean(yTemp(1:JJ))<=yTempKeep
                            yKeep=[yKeep mean(yTemp(1:JJ))];
                            yTempKeep = yKeep(end);
                            xKeep= [xKeep combinedPeakLocs_X{II}];
                            break
                        end
                        
                    end
                end
            end
            
            % now construct a matrix to use for selecting the subset of periods,
            % offsets, and alphas
            selectionMatrix = zeros(length(fillVec),length(ratioVec));
            for II=1:length(xKeep)
                indRow = find(yKeep(II)==fillVec);
                indCol = find(xKeep(II)==ratioVec);
                selectionMatrix(indRow,indCol)=1;
            end
            
            % get subset selections
            obj.extract.periodKeepVec = periodMatrix(find(selectionMatrix));
            obj.extract.angleKeepVec = angleMatrix(find(selectionMatrix));
            obj.extract.scatteringStrengthKeepVec = scatteringStrengthMatrix(find(selectionMatrix));
            obj.extract.fillKeepVec = yKeep;
            obj.extract.ratioKeepVec = xKeep;
            obj.extract.offsetKeepVec = offsetMatrix(find(selectionMatrix));
        end
        
        function obj = synthesize(obj, x, mfd)
            % this method synthesizes a grating given the desired function
            % f
            xo = mfd; 
            gausFunc = @(x,mfd) sqrt(1/(mfd/2))*(2/pi)^(1/4)*exp(-((x-xo)./(mfd/2)).^2);
            g = gausFunc(x,mfd);
            dx = x(2)-x(1); 
            for II=1:length(x)
                alpha(II) = g(II)^2/(1-sum(g(1:II).^2)*dx);                  
            end
            obj.synth.alpha = alpha; 
            figure; hold on
            plot(x,g.^2./max(g.^2))
            plot(x,alpha)
            xlabel('position (um)')
            
            
            
            xPos = 0;
            counter = 1;
            while xPos<x(end)
               P{counter} = obj.returnUnitCell(1000*alpha(idxNearVal(x,xPos)));
               x_position(counter) = xPos; 
               xPos = xPos+P{counter}.period;
               counter = counter+1;
            end           
            obj.synth.P = P; 
            obj.synth.x_position = x_position;             
            
        end
        
        function params = returnUnitCell(obj,alpha)
            % this method returns the parameters of a single unit cell
            % given the desired scattering strength
            ind = idxNearVal(obj.extract.scatteringStrengthKeepVec,alpha);
            
            params.period = obj.extract.periodKeepVec(ind);
            params.angle = obj.extract.angleKeepVec(ind);
            params.fill = obj.extract.fillKeepVec(ind);
            params.ratio = obj.extract.ratioKeepVec(ind);
            params.offset = obj.extract.offsetKeepVec(ind);   
            params.alpha = obj.extract.scatteringStrengthKeepVec(ind);
        end
        
        function obj = generateDimensions(obj)
            % convert obj.synth.P to obj.dims.P real dimensions
            for II=1:length(obj.synth.P)
                period = obj.synth.P{II}.period;
                fill = obj.synth.P{II}.fill;
                ratio = obj.synth.P{II}.ratio;
                offset = obj.synth.P{II}.offset;
                if II==1
                    initP = obj.synth.P{II};
                    dims{II}.period = period;
                    dims{II}.offset = period*offset;
                    dims{II}.cSiWidth = period*fill;
                    dims{II}.pSiWidth = period*fill*ratio;
                    initBool = 1;
                    cntr = 2;
                    % manually duplicate first period
                    dims{cntr} = dims{1};
                    cntr = 3; 
                elseif initBool
                    if ~isequal(obj.synth.P{II},initP)                        
                        dims{cntr}.period = period;
                        dims{cntr}.offset = period*offset;
                        dims{cntr}.cSiWidth = period*fill;
                        dims{cntr}.pSiWidth = period*fill*ratio;
                        cntr = cntr+1;
                        initBool = 0;
                    end
                else
                    dims{cntr}.period = period;
                    dims{cntr}.offset = period*offset;
                    dims{cntr}.cSiWidth = period*fill;
                    dims{cntr}.pSiWidth = period*fill*ratio;
                    cntr = cntr+1;
                end
            end
            
            %             for II=1:length(obj.synth.P)
            %                 period = obj.synth.P{II}.period;
            %                 fill = obj.synth.P{II}.fill;
            %                 ratio = obj.synth.P{II}.ratio;
            %                 offset = obj.synth.P{II}.offset;
            %                 dims{II}.period = period;
            %                 dims{II}.offset = period*offset;
            %                 dims{II}.cSiWidth = period*fill;
            %                 dims{II}.pSiWidth = period*fill*ratio;
            %             end
            
            obj.dims = dims;
        end
        
        function obj = plotParsedMatrices(obj)
            % directivity
            figure;
            imagesc(obj.dataObj.ratioVector, obj.dataObj.fillVector,log10(abs(1./obj.parse.directivityMatrix)))
            set(gca,'ydir','normal'); colorbar; caxis([0 3])
            xlabel('ratio'); ylabel('fill')
            title('directivity')
            
            % period
            figure;
            imagesc(obj.dataObj.ratioVector, obj.dataObj.fillVector,obj.parse.periodMatrix)
            set(gca,'ydir','normal'); colorbar; 
            xlabel('ratio'); ylabel('fill')
            title('period')
            
            % angle
            figure;
            imagesc(obj.dataObj.ratioVector, obj.dataObj.fillVector,obj.parse.angleMatrix)
            set(gca,'ydir','normal'); colorbar; 
            xlabel('ratio'); ylabel('fill')
            title('angle')
            
            % offset
            figure;
            imagesc(obj.dataObj.ratioVector, obj.dataObj.fillVector,obj.parse.offsetMatrix)
            set(gca,'ydir','normal'); colorbar; 
            xlabel('ratio'); ylabel('fill')
            title('offset')
            
            % scattering strength (alpha)
            figure;
            imagesc(obj.dataObj.ratioVector, obj.dataObj.fillVector,obj.parse.scatteringStrengthMatrix)
            set(gca,'ydir','normal'); colorbar; 
            xlabel('ratio'); ylabel('fill')
            title('scattering strength')
        end
        
    end
    
end

