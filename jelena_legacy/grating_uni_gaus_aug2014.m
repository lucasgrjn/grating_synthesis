function gratingTest()
tic
% close all; clear all
design = 'unidirDown1280EOS24Gaussian'; %name of design you want to simulate
dirprefix = strcat(design,'unidirDown1280EOS24Gaussian'); %name of directory where sim files will be stored
dx = 0.005; %x discretization
dy = 0.005; %y discretization
pmlx = 10*dx;
pmly = 10*dy;
P.polarization = 'TE';
P.timeType = 'pulse';
P.numFrames = 40; %number of frames that make up the gif
P.switches.plotSnapshots = 'yes';
P.switches.savePlots = 'no'; 
P.switches.plotData = 'yes';
P.switches.deleteFrames = 'no';  %changed from no to yes
P.switches.justLoadData = 'no'; %changed from no to yes % seems that this does nothing for now
P.switches.addImperfections  = 'no'; % TEM imperfections
P.switches.debug = 'no'; %use this to see if sim is setup correctly before fully running

% for design A
P.l0 = 1.199;
% % for design B
% P.l0 = 1.206;
%P.wavelengthSpectrum = [P.l0-0.02 P.l0+0.02 0.002];    % coarse tune
% P.wavelengthSpectrum = [P.l0-0.01 P.l0+0.01 0.001];    % fine tune
P.wavelengthSpectrum = [P.l0-0.105 P.l0+0.105 0.001];
P.wavelengthVec = P.wavelengthSpectrum(1):P.wavelengthSpectrum(3):P.wavelengthSpectrum(2);
P.wavelengthIndex = find( P.wavelengthVec == P.l0 );
xmax = 12.0; %x domain
ymax = 2.0; %y domain
P.propLength = xmax; %propagation length in um. simulation uses this to determine how long to run sim
if strcmp(P.switches.debug,'yes')
    dirprefix=strcat('DEBUG_',dirprefix);
end
name = strcat(dirprefix);
fdtdInputs = {'numD',2,'userName','Jelena','polarization',P.polarization,'domain', [xmax ymax],...
        'pml', [pmlx pmly],'discretization', [dx dy],'cfl', 1.1,'propLength',P.propLength, ...
        'pulseWidth',50,'pulseDelay',3,'ngGuess', 5,'wavelengthSpectrum',P.wavelengthSpectrum, ...
        'debug', P.switches.debug, ...
        'numMovieFrames', 50,'deleteMovieFrames', P.switches.deleteFrames,'numCores', 7, ...
        'colorScheme', 'redWhiteBlue','namedd', ['dd_',name,'.txt'],'nameDir', name};
gratingSimInputs = {'P',P};
Qgrat = gratingSim_uni_gaus_aug2014(fdtdInputs,gratingSimInputs{:}); %initialize gratingSim object
OPTS.psiAbsorption = 0; %set absorption in dB/cm
OPTS.addImperfections = P.switches.addImperfections; %include TEM imperfections
if strcmp(P.switches.debug,'yes')
    OPTS.OUTPLOTS = [1 zeros(1,10) 0 0 0];
    OPTS.debug = 'yes';
else
    %OPTS.OUTPLOTS =  [1 zeros(1,2) 0 1 0 0 0 0 0 0 0 0 0]; %these choose which plots are generated from getParams
    % OPTS.OUTPLOTS = [ones(1,13) 0];%[1 zeros(1,10) 1 0 0]; %these choose which plots are generated from getParams
    OPTS.OUTPLOTS = [1 1 1 1 1 0 0 0 0 0 0 0 1 0];
    %OPTS.OUTPLOTS = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    OPTS.FSAVEFIGS = 1; OPTS.FIGFILEEXT = 'fig';
end
OPTS.switches.justLoadData = 'no';
Qgrat = Qgrat.runGratingSimulation('process','IBM45','whichDesign',design,'OPTS',OPTS); %run grating simulation
toc