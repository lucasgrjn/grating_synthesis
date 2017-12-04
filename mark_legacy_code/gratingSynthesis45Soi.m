% this script runs the necessary simulation loops to generate the
% structures needed to synthesize a grating radiation pattern

% define simulation variables
d = 0.01; % discretization
dx = d; dy = d;
lambda0 = 1.28;
c0 = 299792458;
mu0 = 4*pi*10^-7;
BC = 0;
numcells = 15; % number of cells to visualize
modes = 20;

% set save directories
saveDir = 'sim/devDir/';
mkdir(saveDir);

% define structure variables
h_csi = .08; %roundToGrid(f_dims_12soi('bodySi','t'),d); %c-Si
h_psi = .08; %roundToGrid(f_dims_12soi('polySi','t'),d); %p-Si
h_sio_bot = 0.12; %roundToGrid(f_dims_12soi('STI','t'),d); %buried oxide
h_sio_top = 0.25; %backend oxide
h_sin = 0.08; %roundToGrid(f_dims_12soi('SiN','t'),d); %silicon nitride liner
h_air = 0.25;
h_pml = 0.5; %pml thickness
h_vec = [h_pml, h_air, h_sio_bot, h_csi, h_psi, h_sin, h_sio_top, h_pml]; %vector of layer thicknesses
h_tot = sum(h_vec);
%refractive index
n_sio_bot = 1.44; %index_IBM12SOI45nm_fits(lambda0,'STI'); %buried oxide
n_sio_top = 1.44; %index_IBM12SOI45nm_fits(lambda0,'PSG'); %backend oxide
n_sin = 1.44; %index_IBM12SOI45nm_fits(lambda0,'SiN'); %silicon nitride liner
n_csi = 3.5; %index_IBM12SOI45nm_fits(lambda0,'bodySi'); %crystalline silicon
n_psi = 3.5; %index_IBM12SOI45nm_fits(lambda0,'polySi'); %polysilicon
n_air = 1.0;

% pml options
PML_options = [1 h_pml 0.6 2];

% fill and layer ratio vectors
fillVec = 0.4:0.02:1.0;
fillVec = 0.4;      % i changed this
ratioVec = 0.1:0.02:1.2;
ratioVec = 0.1;     % i changed this



% define indexes and structures
indexes =   struct(  'background', n_sio_top,...
    'structures', [n_air n_sio_bot n_csi n_sin n_sin n_psi]);



OPTS =      struct( 'numModes', modes, 'BC', BC, 'PML_options', PML_options,...
    'numcells',numcells);
constants = struct( 'omega0', 2*pi/lambda0*2*c0/1e-6, 'mu0', 4*pi*10^-7);

%% generate period matrix
periodVec = 0.46:0.02:0.74;
periodVec = 0.46;       % i changed this
h = waitbar(0,'period search');
numLoops = length(fillVec)*length(ratioVec);
loopCount = 0; 
for II=1:length(fillVec)
    for JJ=1:length(ratioVec)                
        tstart = tic; 
        for KK=1:length(periodVec)
            period = periodVec(KK);
            ratio = ratioVec(JJ);            
            fill = fillVec(II);
            if fill*ratio>1
                ratio = 1/fill;
            end
            offset = 0.0; %set to zero to find period for theta_opt
            
            domain = [period h_tot]; %entire simulation domain including pml
            
            % h_vec = [h_pml, h_air, h_sio_bot, h_csi, h_psi, h_sin, h_sio_top, h_pml]
            refpoints = [0 0; ... % air
                0 sum(h_vec(1:2)); ... % box
                period-fill*period sum(h_vec(1:3)); ... % c-si
                0 sum(h_vec(1:4)); ... % SiN
                period-period*(offset+ratio*fill)-h_sin sum(h_vec(1:4)); ... % SiN
                period-period*(offset+ratio*fill) sum(h_vec(1:4)); ... % p-si
                ];
            
            dimensions = [period h_air+h_pml; period h_sio_bot; period*fill h_csi; period h_sin; ratio*period*fill+2*h_sin h_sin+h_psi; ratio*period*fill h_psi];
            
            
            dims =      struct(  'dims', dimensions,...
                'refpoints', refpoints,...
                'period', period, ...
                'wgRegion', [sum(h_vec(1:2)), sum(h_vec(1:6))] ...
                );
            
            simInputs = {'lambda0', lambda0, 'dims', dims, 'indexes', indexes, 'domain', domain, 'dx', dx, 'dy', dy, 'OPTS', OPTS, 'constants', constants};
            
            gratObj = c_twoLevelGrating_sim(simInputs{:});
            gratObj = gratObj.twoLevelBuilder;
            gratObj = gratObj.runSimulation;
            maxAngle(KK) = gratObj.max_angle;
        end
        periodMatrix(II,JJ) = periodVec(idxNearVal(maxAngle,20));
        loopCount = loopCount+1;
        tstop = toc(tstart);
        waitbar(loopCount/numLoops,h,sprintf('period search, %d seconds left', round(tstop)*(numLoops-loopCount)));
    end
end
close(h)
save([saveDir 'periodMatrix'], 'periodMatrix')

%% find optimal offset to maximize directivity
h = waitbar(0,'offset search');
numLoops = length(fillVec)*length(ratioVec);
loopCount = 0; 
directivityMatrix = zeros(length(fillVec),length(ratioVec));
offsetMatrix = directivityMatrix; 
angleMatrix = directivityMatrix; 
for II=1:length(fillVec)
    for JJ=1:length(ratioVec)
        period = periodMatrix(II,JJ);
        offsetVec = 0:0.015:period;
        tstart=tic;
        directivity = zeros(length(offsetVec));
        maxAngle = zeros(length(offsetVec));
        parfor KK=1:length(offsetVec)            
            ratio = ratioVec(JJ);            
            fill = fillVec(II);
            if fill*ratio>1
                ratio = 1/fill
            end
            offset = offsetVec(KK);
            
            domain = [period h_tot]; %entire simulation domain including pml
            
            % h_vec = [h_pml, h_air, h_sio_bot, h_csi, h_psi, h_sin, h_sio_top, h_pml]
            refpoints = [0 0; ... % air
                0 sum(h_vec(1:2)); ... % box
                period-fill*period sum(h_vec(1:3)); ... % c-si
                0 sum(h_vec(1:4)); ... % SiN
                period-period*(offset+ratio*fill)-h_sin sum(h_vec(1:4)); ... % SiN
                period-period*(offset+ratio*fill) sum(h_vec(1:4)); ... % p-si
                ];
            
            dimensions = [period h_air+h_pml; period h_sio_bot; period*fill h_csi; period h_sin; ratio*period*fill+2*h_sin h_sin+h_psi; ratio*period*fill h_psi];
            
            
            dims =      struct(  'dims', dimensions,...
                'refpoints', refpoints,...
                'period', period, ...
                'wgRegion', [sum(h_vec(1:2)), sum(h_vec(1:6))] ...
                );
            
            simInputs = {'lambda0', lambda0, 'dims', dims, 'indexes', indexes, 'domain', domain, 'dx', dx, 'dy', dy, 'OPTS', OPTS, 'constants', constants};
            
            gratObj = c_twoLevelGrating_sim(simInputs{:});
            gratObj = gratObj.twoLevelBuilder;
            gratObj = gratObj.runSimulation;
            directivity(KK) = gratObj.directivity;
            maxAngle(KK) = gratObj.max_angle; 
        end
        [directivityMatrix(II,JJ) idxMinDirectivity] = min(directivity);
        offsetMatrix(II,JJ) = offsetVec(idxMinDirectivity); 
        angleMatrix(II,JJ) = maxAngle(idxMinDirectivity); 
        loopCount = loopCount+1;
        tstop = toc(tstart);
        waitbar(loopCount/numLoops,h,sprintf('offset search, %d seconds left', round(tstop)*(numLoops-loopCount)));
    end
end
save([saveDir 'directivityMatrix'], 'directivityMatrix')
save([saveDir 'angleMatrix'], 'angleMatrix')
save([saveDir 'offsetMatrix'], 'offsetMatrix')
close(h)






