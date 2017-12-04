classdef c_gratObj45Soi
    %C_GRATOBJ45SOI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fillVector
        ratioVector
        periodVector
        offsetVector
        periodMatrix %note: these are really 4 dimensional arrays. need to change naming from Matrix to Array but will break some things
        offsetMatrix
        scatteringStrengthMatrix
        directivityMatrix
        angleMatrix
        lambda0
        P
        inputs
    end
    
    methods
        
        function obj = c_gratObj45Soi(varargin)
            domainFields = {'lambda0','d','fillVector','ratioVector','periodVector','offsetVector','optimal_angle'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case domainFields
                        eval(['obj.inputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            obj.fillVector = obj.inputs.fillVector; 
            obj.ratioVector = obj.inputs.ratioVector; 
            obj.periodVector = obj.inputs.periodVector; 
            obj.offsetVector = obj.inputs.offsetVector; 
        end
        
        function obj = setParams(obj)
            % extract some inputs
            obj.P.d = obj.inputs.d;
            obj.P.lambda0 = obj.inputs.lambda0;
            
            % set layer thicknesses
            obj.P.thickness.csi = roundToGrid(f_dims_12soi('bodySi','t'),obj.P.d); %crystalline Si
            obj.P.thickness.psi = roundToGrid(f_dims_12soi('polySi','t'),obj.P.d); %p-Si
            obj.P.thickness.sio_bot = roundToGrid(f_dims_12soi('STI','t'),obj.P.d); %buried oxide
            obj.P.thickness.sio_top = 0.25; %backend oxide
            obj.P.thickness.sin = roundToGrid(f_dims_12soi('SiN','t'),obj.P.d); %silicon nitride liner
            obj.P.thickness.air = 0.25;
            obj.P.thickness.pml = 0.5; %pml thickness
            obj.P.thickness.vec = [obj.P.thickness.pml, obj.P.thickness.air, obj.P.thickness.sio_bot, obj.P.thickness.csi, obj.P.thickness.psi, obj.P.thickness.sin, obj.P.thickness.sio_top, obj.P.thickness.pml]; %vector of layer thicknesses
            obj.P.thickness.sum = sum(obj.P.thickness.vec); % total thickness
            
            % set refractive indexes
            obj.P.index.sio_bot = index_IBM12SOI45nm_fits(obj.P.lambda0,'STI'); %buried oxide
            obj.P.index.sio_top = index_IBM12SOI45nm_fits(obj.P.lambda0,'PSG'); %backend oxide
            obj.P.index.sin = index_IBM12SOI45nm_fits(obj.P.lambda0,'SiN'); %silicon nitride liner
            obj.P.index.csi = index_IBM12SOI45nm_fits(obj.P.lambda0,'bodySi'); %crystalline silicon
            obj.P.index.psi = index_IBM12SOI45nm_fits(obj.P.lambda0,'polySi'); %polysilicon
            obj.P.index.air = 1.0;
            
            % pml options
            PML_options = [1 obj.P.thickness.pml 0.6 2];
            
            % set index structures
            obj.P.index_structures = struct(  'background', obj.P.index.sio_top,...
                'structures', [obj.P.index.air obj.P.index.sio_bot obj.P.index.csi obj.P.index.sin obj.P.index.sin obj.P.index.psi]);
            
            %set options
            c0 = 299792458;
            mu0 = 4*pi*10^-7;
            BC = 0;
            numcells = 15; % number of cells to visualize
            modes = 20;
            obj.P.OPTS = struct( 'numModes', modes, 'BC', BC, 'PML_options', PML_options,...
                'numcells',numcells);
            obj.P.constants = struct('omega0', 2*pi/obj.P.lambda0*2*c0/1e-6, 'mu0', 4*pi*10^-7);
            
        end
        
        function obj = runParameterSweep(obj)
            % this method runs a full parameter sweep
            
            % extract some variables from the object
            fillVec = obj.fillVector;
            ratioVec = obj.ratioVector;
            periodVec = obj.periodVector;
            offsetVec = obj.offsetVector;
            
            thickness_vec = obj.P.thickness.vec;
            
            % run loops
            numLoops = length(fillVec)*length(offsetVec)*length(ratioVec)*length(periodVec);
            h = waitbar(0, 'Loops running. God help us all.');
            loopCount = 0;
            
            periodMatrix = zeros(length(fillVec),length(ratioVec),length(offsetVec),length(periodVec));
            offsetMatrix = periodMatrix;
            scatteringStrengthMatrix = periodMatrix;
            directivityMatrix= periodMatrix;
            angleMatrix = periodMatrix;
            
            for II=1:length(fillVec)
                fill = fillVec(II);
                for JJ=1:length(ratioVec)
                    ratio = ratioVec(JJ);
                    if fill*ratio>1 % make sure poly isn't larger than period
                        ratio = 1/fill;
                    end
                    
                    for KK=1:length(offsetVec)
                        offset = offsetVec(KK);
                        tstart = tic;
                        parfor LL=1:length(periodVec)
                            period = periodVec(LL);
                            domain = [period obj.P.thickness.sum];
                            
                            refpoints = [0 0; ... % air
                                0 sum(thickness_vec(1:2)); ... % box
                                period-fill*period sum(thickness_vec(1:3)); ... % c-si
                                0 sum(thickness_vec(1:4)); ... % SiN
                                period-period*(offset+ratio*fill)-obj.P.thickness.sin sum(thickness_vec(1:4)); ... % SiN
                                period-period*(offset+ratio*fill) sum(thickness_vec(1:4)); ... % p-si
                                ];
                            
                            dimensions = [period obj.P.thickness.air+obj.P.thickness.pml; ...
                                period obj.P.thickness.sio_bot; ...
                                period*fill obj.P.thickness.csi; ...
                                period obj.P.thickness.sin; ...
                                ratio*period*fill+2*obj.P.thickness.sin obj.P.thickness.sin+obj.P.thickness.psi; ...
                                ratio*period*fill obj.P.thickness.psi ...
                                ];
                            
                            dims =  struct(     'dims', dimensions, ...
                                'refpoints', refpoints, ...
                                'period', period, ...
                                'wgRegion', [sum(thickness_vec(1:2)), sum(thickness_vec(1:6))] ...
                                );
                            
                            simInputs = {'lambda0',obj.P.lambda0, 'dims',dims, 'indexes', obj.P.index_structures, 'domain', domain, 'dx', obj.P.d, 'dy', obj.P.d, 'OPTS', obj.P.OPTS, 'constants', obj.P.constants};
                            
                            tempObj = c_twoLevelGrating_sim(simInputs{:});

                            tempObj = tempObj.twoLevelBuilder;

                            tempObj = tempObj.runSimulation;
                            
                            periodMatrix(II,JJ,KK,LL) = period;
                            offsetMatrix(II,JJ,KK,LL) = offset;
                            scatteringStrengthMatrix(II,JJ,KK,LL) = tempObj.alpha;
                            directivityMatrix(II,JJ,KK,LL) = tempObj.directivity;
                            angleMatrix(II,JJ,KK,LL) = tempObj.max_angle;
                            
                        end
                        tstop = toc(tstart);
                        loopCount=loopCount+1;
                        waitbar(loopCount/numLoops*length(periodVec), h, sprintf('Loops running, %d seconds left', round(tstop*(numLoops/length(periodVec)-loopCount))));
                    end
                end
            end
            obj.periodMatrix = periodMatrix;
            obj.offsetMatrix = offsetMatrix;
            obj.scatteringStrengthMatrix = scatteringStrengthMatrix;
            obj.directivityMatrix = directivityMatrix;
            obj.angleMatrix = angleMatrix;
            close(h);
        end              
        
    end
    
end

