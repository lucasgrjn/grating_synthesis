% test the c_synthesizeGrating class

fileName = 'sim/gratObjFineSweep.mat';
optimalAngle = 20.0; 

inputs = {'fileName',fileName,'optimalAngle',optimalAngle};

% initialize
synthObj = c_synthesizeGrating(inputs{:});

% load
synthObj = synthObj.loadData;
% parse
synthObj = synthObj.parseData;
% extract contours
synthObj = synthObj.extractContours; 
%% synthesize field
mfd = 10.0;
x = 0:0.1:3*mfd; 
synthObj = synthObj.synthesize(x,mfd);
synthObj = synthObj.generateDimensions;

for II=1:length(synthObj.synth.P)
    offsets(II) = synthObj.synth.P{II}.offset;    
    alphas(II) = synthObj.synth.P{II}.alpha;  
    periods(II) = synthObj.synth.P{II}.period; 
    fills(II) = synthObj.synth.P{II}.fill; 
    ratios(II) = synthObj.synth.P{II}.ratio;
    angles(II) = synthObj.synth.P{II}.angle;
end
xPos = synthObj.synth.x_position; 
%plot(xPos,alphas./1e3,'-o')
%legend({'Gaussian profile' 'scattering strength' 'synthesized'})
%figure; hold on
%plot(xPos,alphas./1e3,'-o')
%plot(xPos,offsets,'-o')
%plot(xPos,periods,'-o')
%plot(xPos,fills,'-o')
%plot(xPos,ratios,'-o')
%legend({'alpha' 'offset' 'period' 'fill' 'ratio'})
%xlabel('position (um)')
%figure; plot(xPos,angles,'-o'); xlabel('position (um)'); ylabel('angle (degrees)')
%synthObj.plotParsedMatrices;

save('synthObj.mat','synthObj');