% NOTE: this is a development stage script! beware of bugs!
% it's hard coded to find the downward designs!

% load the sweep data

tmp = load('sim/gratObjFineSweep.mat');
gratObj = tmp.gratObj;

% get some variables
fillVec = gratObj.fillVector;
ratioVec = gratObj.ratioVector;
offsetVec = gratObj.offsetVector;
periodVec = gratObj.periodVector;

% optimal angle
optimalAngle = 20.0; %degrees

% G(fill, ratio, off, period)

% for each fill, ratio point, find the period that gives closest to 20
% degrees output radiation

for II=1:length(fillVec)
    for JJ=1:length(ratioVec)
        % for each offset (row), find period that gives angle closest to
        % optimal angle
        tempAngleMatrix = squeeze(gratObj.angleMatrix(II,JJ,:,:));
        [minErrorVec, minAngIdx] = min(abs(tempAngleMatrix-optimalAngle),[],2);
        dMatrix = squeeze(gratObj.directivityMatrix(II,JJ,:,:));
        
        for KK=1:length(minAngIdx)
            dVec(KK) = dMatrix(KK,minAngIdx(KK));            
        end
        
        [minD, minDIdx] = min(abs(dVec));
        
        % keep these variables      
        try
            offsetMatrix(II,JJ) = offsetVec(minDIdx);
        catch
            disp([II JJ])
        end
        periodMatrix(II,JJ) = periodVec(minAngIdx(minDIdx));
        angleMatrix(II,JJ) = tempAngleMatrix(minDIdx,minAngIdx(minDIdx));
        directivityMatrix(II,JJ) = minD;
        scatteringStrengthMatrix(II,JJ) = gratObj.scatteringStrengthMatrix(II,JJ,minDIdx,minAngIdx(minDIdx));
        
    end
end

%% figure out directivity contour extraction
%peakLocs = nan(length(ratioVec),2);
% for II=1:length(ratioVec)
%     dVecLine = 1./directivityMatrix(:,II);
%     [~,locs]=findpeaks(dVecLine,fillVec,'NPeaks',2,'MinPeakHeight',10,'MinPeakProminence',10);
%     if length(locs)>1
%         peakLocs(II,:)=locs;
%     else
%         peakLocs(II,1)=locs;
%     end    
% end
% scan over fillVec (y)
peakLocs = nan(length(fillVec),4);
logDirectivityDown = log10(abs(1./directivityMatrix)); 
myFilter = fspecial('gaussian',[3 3], 1.0);
filteredLogDirectivityDown = imfilter(logDirectivityDown, myFilter, 'replicate');

for II=1:length(fillVec)
    %dVecLine = log10(abs(1./directivityMatrix(II,:)));
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
figure; hold on
for II=1:size(peakLocs,2)
    plot(peakLocs(:,II),fillVec,'ro')
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
            try
                peakLocY_1{JJ,II}=tempPeakLocY(idx,II);
            catch
                disp('here')
            end
        end
    end
end


%% scan over ratioVec (x)
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
%figure; hold on
for II=1:size(peakLocs,2)
    plot(ratioVec,peakLocs(:,II),'go')
end

% parse through the data to extract points to fit 
peakLocX_2 = cell(length(ratioVec),size(peakLocs,2)); %cell array
peakLocY_2 = cell(length(ratioVec),size(peakLocs,2)); %cell array

for II=1:size(peakLocs,2)
%     tempVec = peakLocs(II,:);
%     [sortedVec sortedIdx] = sort(tempVec);
%     tempPeakLocX_2(:,II) = sortedVec; 
%     tempPeakLocY_2(:,II) = fillVec(sortedIdx);
%     for JJ=1:size(peakLocs,1)
%         idx = find(tempPeakLocX_2(:,II)==ratioVec(JJ));
%         peakLocX_2{JJ,II}=ratioVec(JJ);
%         if isempty(idx)
%             peakLocY_2{JJ,II}=nan;
%         else
%             try
%                 peakLocY_2{JJ,II}=tempPeakLocY_2(idx,II);
%             catch
%                 disp('here')
%             end
%         end
%     end
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

%%
% now construct a matrix to use for selecting the subset of periods,
% offsets, and alphas
selectionMatrix = zeros(length(fillVec),length(ratioVec));
for II=1:length(xKeep)
    indRow = find(yKeep(II)==fillVec);
    indCol = find(xKeep(II)==ratioVec);
    selectionMatrix(indRow,indCol)=1;
end

imagesc(selectionMatrix); set(gca,'ydir','normal'); 

%% get subset selections
periodKeepVec = periodMatrix(find(selectionMatrix));
scatteringStrengthKeepVec = scatteringStrengthMatrix(find(selectionMatrix));
fillKeepVec = yKeep; 
ratioKeepVec = xKeep;
offsetKeepVec = offsetMatrix(find(selectionMatrix));

plot(ratioKeepVec,scatteringStrengthKeepVec,'o');



%% fit and plot the inverted design curve
fit = polyfit(xKeep,yKeep,5);
figure;
plot(xKeep,yKeep,'o',xKeep,polyval(fit,xKeep))





%%
%figure; imagesc(ratioVec,fillVec,scatteringStrengthMatrix); set(gca,'ydir','normal'); colorbar; title('scattering strength'); xlabel('ratio'); ylabel('fill')
%figure; imagesc(ratioVec,fillVec,log10(abs(1./directivityMatrix))); set(gca,'ydir','normal'); colorbar; title('downward directivity'); xlabel('ratio'); ylabel('fill')
%figure; imagesc(ratioVec,fillVec,angleMatrix); set(gca,'ydir','normal'); colorbar; title('angle'); xlabel('ratio'); ylabel('fill')
%figure; imagesc(ratioVec,fillVec,periodMatrix); set(gca,'ydir','normal'); colorbar; title('period'); xlabel('ratio'); ylabel('fill')
%figure; imagesc(ratioVec,fillVec,offsetMatrix); set(gca,'ydir','normal'); colorbar; title('offset'); xlabel('ratio'); ylabel('fill')

%plot fit line on top
figure; hold on
imagesc(ratioVec,fillVec,log10(abs(1./directivityMatrix))); set(gca,'ydir','normal'); colorbar; title('downward directivity'); xlabel('ratio'); ylabel('fill')
plot(xKeep,polyval(fit,xKeep))


logDirectivityDown = log10(abs(1./directivityMatrix)); 
myFilter = fspecial('gaussian',[3 3], 1.0);
filteredLogDirectivityDown = imfilter(logDirectivityDown, myFilter, 'replicate');

figure; hold on
imagesc(ratioVec,fillVec,filteredLogDirectivityDown); set(gca,'ydir','normal'); colorbar; title('downward directivity'); xlabel('ratio'); ylabel('fill')
f_placeContoursOnPlot(ratioVec,fillVec,filteredLogDirectivityDown,[1:3]);


%localMaxMat = localmax(filteredLogDirectivityDown);
%figure; imagesc(ratioVec,fillVec,localMaxMat); set(gca,'ydir','normal'); colorbar; title('downward directivity'); xlabel('ratio'); ylabel('fill')


