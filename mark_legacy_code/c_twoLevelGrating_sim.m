classdef c_twoLevelGrating_sim
    % encapsulates a 2 level grating simulation
    
    properties
        OPTS
        dims
        indexes
        N
        domain
        dx
        dy
        inputs
        constants
        max_angle
        directivity
        alpha
        k
        E_z
    end
    
    methods
        function obj = c_twoLevelGrating_sim(varargin)
            domainFields = {'dims', 'OPTS', 'indexes', 'domain', 'dx', 'dy', 'lambda0','constants'};
            for k = 1:2:size(varargin,2)
                varName = varargin{k};
                switch(varName)
                    case domainFields
                        eval(['obj.inputs.' varName ' = varargin{k+1};']);
                    otherwise
                        error(['Trying to assign value to invalid parameter "' varargin{k} '".  Check grating simulation instructions.']);
                end
            end
            % check some inputs
            obj.OPTS = obj.inputs.OPTS;
            obj.dims = obj.inputs.dims;
            obj.indexes = obj.inputs.indexes;
            obj.domain = obj.inputs.domain;
            obj.dx = obj.inputs.dx;
            obj.dy = obj.inputs.dy;
            obj.constants = obj.inputs.constants;
        end
        
        function obj = twoLevelBuilder(obj)
            domain = obj.domain; dx = obj.dx; dy = obj.dy; refpoint = obj.dims.refpoints;
            dims = obj.dims.dims; % kinda bad naming scheme...
            nlyrs = obj.indexes.structures; nbackground = obj.indexes.background;
            
            N = nbackground*ones(round(domain(1)/dx)+1, round(domain(2)/dy)+1).';
            
            % discretize inputs %NOTE: feb 21, 2016 changed from midpoints
            % to bottom left corner and calling it "refpoints"
            refpoints_d = [round(refpoint(:,1)./dx), round(refpoint(:,2)./dy)];
            dims_d = [round(dims(:,1)./dx), round(dims(:,2)./dy)];
            
            % place structures
            for II=1:length(nlyrs)
%                 startIndX = max(0,floor(refpoints_d(II,1)))+1;
%                 endIndX = min(size(N,2),floor(refpoints_d(II,1)+dims_d(II,1))+1);
%                 startIndY = max(0,round(refpoints_d(II,2)))+1; %currently not treating wrapping y structures
%                 endIndY = min(size(N,1),round(refpoints_d(II,2)+dims_d(II,2)))+1;
%                 N(startIndY:endIndY,startIndX:endIndX)=nlyrs(II);
                %below is legacy from a plot I was generating for my thesis
                startIndX = floor(refpoints_d(II,1))+1;
                endIndX = floor(refpoints_d(II,1)+dims_d(II,1))+1;
                w_pixels = endIndX-startIndX; % width in pixels
                if startIndX<=0 && w_pixels<size(N,2) %split into two structures
                    SPLITBOOL=1;
                    startIndXLeft = 1; % left structure
                    endIndXLeft = endIndX+1; % left structure
                    startIndXRight = size(N,2)+startIndX; % right structure
                    endIndXRight = size(N,2); % right structure
                    %sum([endIndXLeft endIndXRight] - [startIndXLeft startIndXRight])
                    %w_pixels
                elseif endIndX>size(N,2) && w_pixels<size(N,2)
                    SPLITBOOL=1;
                    startIndXLeft = 1;
                    endIndXLeft = endIndX-size(N,2)+1;
                    startIndXRight = startIndX;
                    endIndXRight = size(N,2);
                    %sum([endIndXLeft endIndXRight] - [startIndXLeft startIndXRight])
                    %w_pixels
                elseif w_pixels>=size(N,2)
                    SPLITBOOL=0;
                    startIndX=1;
                    endIndX=size(N,2);
                else
                    SPLITBOOL=0;
                end
                startIndY = round(refpoints_d(II,2))+1; %currently not treating wrapping y structures
                endIndY = round(refpoints_d(II,2)+dims_d(II,2))+1;
                if SPLITBOOL
                    N(startIndY:endIndY,startIndXLeft:endIndXLeft) = nlyrs(II);
                    N(startIndY:endIndY,startIndXRight:endIndXRight) = nlyrs(II);
                else
                    N(startIndY:endIndY,startIndX:endIndX) = nlyrs(II);
                end
            end
            obj.N = N;
        end
        
        function obj = runSimulation(obj)
            a = obj.dims.period;
            % run solver
            guessk = pi/(2*a);
            [Phi_1D, k] = complexk_mode_solver_2D_PML(obj.N,obj.dx,2*pi/obj.inputs.lambda0,obj.OPTS.numModes,guessk,obj.OPTS.BC,obj.OPTS.PML_options);
            % sort modes based on guided power
            number_x = size(obj.N,2);
            number_y = size(obj.N,1);
            guidepower = zeros(length(k),1);
            Phi = zeros(number_y, number_x);
            for II = 1:length(k)
                for JJ = 1:number_y
                    Phi(JJ,:) = Phi_1D((JJ-1)*number_x+1:JJ*number_x,II);
                end
                bot = round(obj.dims.wgRegion(1)/obj.dy);
                top = round(obj.dims.wgRegion(2)/obj.dy);
                guidepower(II) = sum(sum(abs(Phi(bot:top,:))))/sum(sum(abs(Phi))); %estimate of power contained in waveguide
            end
            sortedmatrix = sortrows([guidepower k Phi_1D.'],-1);
            k = sortedmatrix(1,2);
            % only keep decaying k
            if imag(k) >= 0
                numcells = obj.OPTS.numcells;
                Phi_1D = (sortedmatrix(1,3:end)).';
                for JJ=1:number_y
                    Phi(JJ,:) = Phi_1D((JJ-1)*number_x+1:JJ*number_x);
                end
                E_z = zeros(number_y,number_x*numcells);
                for JJ=0:numcells-1
                    E_z(:,(1+JJ*number_x):(number_x+JJ*number_x)) = Phi;
                end
                for JJ = 1:size(E_z,2)
                    E_z(:,JJ) = E_z(:,JJ)*exp(1i*sortedmatrix(1,2)*obj.dy*JJ);
                end
                % power radiated up
                h_pml_d = round(obj.OPTS.PML_options(2)/obj.dy); %PML_options(2) = h_pml (height of pml)
                y_up = number_y-h_pml_d+1;
                omega0 = obj.constants.omega0; mu0 = obj.constants.mu0;
                H_x_up = 1/(1i*omega0*mu0) .* (E_z(y_up+1,:)-E_z(y_up-1,:))/(2*obj.dx);
                P_rad_up = sum(real(H_x_up(:)'.*E_z(y_up,:)))*obj.dx;
                % power radiated down
                y_dn = h_pml_d;
                H_x_dn = 1/(1i*omega0*mu0) .* (E_z(y_dn+1,:)-E_z(y_dn-1,:))/(2*obj.dx);
                P_rad_dn = -sum(real(H_x_dn(:)'.*E_z(y_dn,:)))*obj.dx;
                % power in (at left edge)
                H_y_in = 1/(-1i*omega0*mu0) .* (E_z(y_dn:y_up,2)-E_z(y_dn:y_up,1))/(obj.dy);
                P_in = -sum(real(conj(H_y_in(:)).*E_z(y_dn:y_up,1)))*obj.dx;
                % calculate alpha efficiency (in 1/um)
                P_rad_dn_unitcell = -sum(real(conj(H_x_dn(1:number_x)).*E_z(y_dn,1:number_x)))*obj.dx;
                obj.alpha = P_rad_dn_unitcell/(P_in*a);
                % directivity
                obj.directivity = P_rad_up/P_rad_dn;
                obj.k = k;
                % Other variables
                k_y_dn = sqrt((2*pi/obj.inputs.lambda0*obj.N(h_pml_d+1,1))^2-k^2);
                obj.max_angle = atand(real(k)/real(k_y_dn));
                obj.E_z = E_z;
            else
                obj.directivity = nan;
                obj.max_angle = nan;
                obj.alpha = nan;
            end
        end
        
        function plotIndex(obj)
            figure
            h=imagesc([1:(size(obj.N,2)-1)].*obj.dx, [1:(size(obj.N,1)-1)]*obj.dy,obj.N);
            set(gca, 'ydir','normal')
            axis image; colormap(h.Parent,flipud(gray)); caxis(h.Parent,[1 max(obj.N(:))]); colorbar
        end
        
        function plotEz(obj)
            imagesc(abs(obj.E_z))
        end
        
    end
    
end

