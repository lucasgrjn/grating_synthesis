% test simulation of 2 level grating
%close all; clear all 

% define structure variables
h_si_bot = 0.080; %c-Si 
h_si_top = 0.080; %p-Si
h_sio_bot = 2.0; 
h_sio_top = 2.0;
h_pml = 0.5; 
%a = 0.4:0.05:0.6;  % period
a = 0.45:0.05:0.6;
%a = 0.54:0.005:0.56; 
a = 0.4:0.02:0.6;
a = 0.5; 
%fill = [0.75:0.005:0.76]; % silicon width/period
fill = 0.75:0.005:0.76;
fill = 0.7:0.05:0.8;
fill = 0.2:0.1:0.8; 
fill = 0.7;
%fill = 0.74:0.005:0.76;
n_sio = 1.5; 
n_si = 3.5;
h_tot = sum([h_pml, h_si_bot, h_si_top, h_sio_bot, h_sio_top, h_pml]);

% define simulation variables
d = 0.01; % discretization 
dx = d; dy = d; % Jelena's code doesn't support different dx,dy!
lambda0 = 1.18;
c0 = 299792458;
mu0 = 4*pi*10^-7;
PML_options = [1 h_pml 0.6 2];
BC = 0; 
numcells = 30; % number of cells to visualize
modes = 20;



% define indexes and structures
indexes =   struct(  'background', 1.5,...
                        'structures', [n_si n_si]); 


                    
OPTS =      struct( 'numModes', modes, 'BC', BC, 'PML_options', PML_options,...
                    'numcells',numcells);
constants = struct( 'omega0', 2*pi/lambda0*2*c0/1e-6, 'mu0', 4*pi*10^-7);                 



%%

% calculate maximum offsets to span
% offsets

numOffsets = 3;

for kk=1:length(a) 
    domain = [a(kk) h_tot];
    for ii=1:length(fill)
        %offset = floor(-a/2/d)*d:d:ceil(a/2/d)*d;
        offset = 0.11
        %offset = 0:d:ceil(a/2/d)*d;
        for jj=1:length(offset)
            
            disp([kk ii jj])
            
%             dims =      struct(  'dims', roundToGrid([a(kk)*fill(ii) h_si_top; a(kk)*fill(ii) h_si_bot],d),...
%                 'midpoints', floor([a(kk)/2-offset(jj)/2 h_tot/2; a(kk)/2+offset(jj)/2 h_tot/2-h_si_bot]./d).*d,...
%                 'period', a(kk), ...
%                 'wgRegion', [h_pml+h_sio_bot, h_pml+h_sio_bot+h_si_bot+h_si_top]);
            dims =      struct(  'dims', roundToGrid([a(kk)*fill(ii) h_si_top; a(kk)*fill(ii) h_si_bot],d),...
                'midpoints', round([0 h_tot/2; offset(jj) h_tot/2-h_si_bot]./d).*d,...
                'period', a(kk), ...
                'wgRegion', [h_pml+h_sio_bot, h_pml+h_sio_bot+h_si_bot+h_si_top]);

                        
            simInputs = {'lambda0', lambda0, 'dims', dims, 'indexes', indexes, 'domain', domain, 'dx', dx, 'dy', dy, 'OPTS', OPTS, 'constants', constants};
            
%             gratObj{KK,II,JJ} = c_twoLevelGrating_sim(simInputs{:});
%             gratObj{KK,II,JJ} = gratObj{KK,II,JJ}.twoLevelBuilder;
%             gratObj{KK,II,JJ} = gratObj{KK,II,JJ}.runSimulation;
%             gratObj{KK,II,JJ}.directivity
            gratObj = c_twoLevelGrating_sim(simInputs{:});
            gratObj = gratObj.twoLevelBuilder;
            gratObj = gratObj.runSimulation;
            gratObj.directivity

            directivity(kk,ii,jj) = gratObj.directivity;
            gratObjStore{kk,ii,jj} = gratObj; 
            
        end
    end
end
%% plot radiated field
fld = gratObj.E_z;
fld = fld/max(fld(:));
x = (1:size(gratObj.N,2)*numcells).*gratObj.dx;
y = (1:size(gratObj.N,1)).*gratObj.dy;

%figure; imagesc(real(fldAngle)); caxis([-1 1])
figure; imagesc(x,y,flipud(real(fld))); colormap(redbluehilight); caxis([-0.5 0.5])
xlabel('x (um)'); ylabel('y (um)'); set(gca,'ydir','normal')
axis image; set(gca,'FontSize',14); set(gcf,'color',[1 1 1])



%% plot directivity
[val,ind] = max(directivity(:));
[r,c,z] = ind2sub(size(directivity),ind);
a(r)
fill(c)
offset(z)

plot(offset,log10(squeeze(directivity(r,c,:))),'o--')


