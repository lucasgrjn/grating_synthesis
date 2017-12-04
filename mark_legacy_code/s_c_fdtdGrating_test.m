% tests the functionality of c_fdtdGrating.m

%% fdtd inputs

close all; clear all
dirprefix = strcat('sim/fdtd/','test'); %name of directory where sim files will be stored

dx = 0.01; %x discretization
dy = 0.01; %y discretization
pmlx = 10*dx;
pmly = 10*dy;

P.polarization = 'TE';
P.propLength = 25; %propagation length in um. simulation uses this to determine how long to run sim
P.pulseWidth = 50; 
P.pulseDelay = 3; 

P.timeType = 'pulse';
P.numMovieFrames = 10; %number of frames that make up the gif
P.switches.plotSnapshots = 'yes';
P.switches.savePlots = 'yes'; 
P.switches.plotData = 'yes';
P.switches.deleteFrames = 'no';  %changed from no to yes
P.switches.justLoadData = 'no'; %changed from no to yes % seems that this does nothing for now
P.switches.debug = 'no'; %use this to see if sim is setup correctly before fully running

P.l0 = 1.280;
P.wavelengthSpectrum = [P.l0-0.1 P.l0+0.1 0.001];
P.wavelengthVec = P.wavelengthSpectrum(1):P.wavelengthSpectrum(3):P.wavelengthSpectrum(2);
P.wavelengthIndex = find( P.wavelengthVec == P.l0 );

if strcmp(P.switches.debug,'yes')
    dirprefix=strcat('DEBUG_',dirprefix);
end    

name = strcat(dirprefix);

%%

% load synthObj
fileName = 'synthObj.mat';
obj = load(fileName);
obj = obj.synthObj; 
P.dims = obj.dims; 

fdtdInputs = {'numD',2,...
        'userName','markwade',...
        'polarization',P.polarization,...
        'domain', [0 0],... % placeholder... gets set in c_fdtdGrating
        'pml', [pmlx pmly], ...
        'discretization', [dx dy], ...
        'cfl', 1.1, ...
        'propLength',P.propLength, ...
        'pulseWidth',P.pulseWidth, ...
        'pulseDelay',P.pulseDelay,...
        'ngGuess', 5, ...%         'backgroundIndex', L.sio_index_1550, ...
        'wavelengthSpectrum',P.wavelengthSpectrum, ...
        'debug', P.switches.debug, ...
        'numMovieFrames', P.numMovieFrames, ...
        'deleteMovieFrames', P.switches.deleteFrames, ...
        'numCores', 6, ...
        'codePath', '/Users/markwade/Documents/Code/fdtd/bin/osx', ...
        'batchfilepath', '/Users/markwade/Documents/Code/fdtd/shell', ...
        'colorScheme', 'redWhiteBlue', ...
        'namedd', ['dd_',name,'.txt'], ...
        'nameDir', name};
    
gratingSimInputs = {'P',P};
Qgrat = c_fdtdGrating(fdtdInputs,gratingSimInputs{:}); %initialize gratingSim object
Qgrat = Qgrat.getParams; 
Qgrat = Qgrat.createGrating; 
Qgrat = Qgrat.runSimulation;
Qgrat = Qgrat.processData('OPTS',struct());

Qgrat.plotDielectrics;
    

    
    
    



