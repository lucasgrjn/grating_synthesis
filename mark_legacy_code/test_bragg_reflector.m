% testing jelena's bloch code

% trying to recreate results from fig 1 in opt lett. paper
% which is a 2 layer quarterwave bragg reflector

% Jelena's FDFD 2D complex-k mode solver
% Version 3
% May 15, 2014
% Includes: Perfect Electric or Magnetic Boundary Conditions, Periodic
% Boundary Conditions, Pefectly Matched Layers (PMLs)
%
% INPUTS:
% n: matrix of refractive indices for the structure
% d: discretization size (nm)
% k0: free-space wavevector (1/nm)
% modes: number of modes to find
% guessk: guess value of complex k
% BC: y boundary condition (PMC=1 or PEC=0) 
% PML_options(1): PML in y direction (yes=1 or no=0)
% PML_options(2): length of PML layer in nm
% PML_options(3): strength of PML in the complex plane
% PML_options(4): PML polynomial order (1, 2, 3...)
%
% OUTPUTS:
% Phi_1D: matrix of all Phi-field solutions (E = e^(jkx)*Phi) in 1D form for desired structure
% k: vector of complex wavevector values calulated

% all units nm

close all; clc; clear;

% -------------------------------------------------------------------------
% PARAMETERS

% indices
n1 = 1;
n2 = 2;

% here I want the period to be a perfect integer number of steps
% pick discretization first
dxz = 1; % nm, same discretization in x and z

m = round(1500/ (4*dxz/(1/n1+1/n2)));

% wavelength
lambda0 = m*dxz*4/(1/n1+1/n2);     % in nm, see notes for why i chose this
k0 = 2*pi/lambda0;

% now i can calculate dimensions of entire space
X = 5;
z1 = lambda0/(4*n1);
z2 = lambda0/(4*n2);
a = z1 + z2; % bloch period
Z = a;

% make x & z space vectors
Z/dxz
nX = round(X/dxz); nZ = round(Z/dxz); % number of data points
x_vec = (-nX/2:nX/2-1).*dxz;
z_vec = (0:nZ-1).*dxz;

% -------------------------------------------------------------------------
% INDEX UNIT CELL

% index refraction
n = n1*ones(nX,nZ);     % dimensions (x, z)
n( :, z_vec >= z1 ) = n2;

% display
figure; imagesc(n); set(gca,'Ydir','normal'); colorbar;
xlabel('z'); ylabel('x');

% -------------------------------------------------------------------------
% RUN MODESOLVER

% inputs
k0 = 1/a;           % testing somethin
% guessk = 5*pi/(8*a);
modes = 1;
BC = 1; % 0 = PEc, 1 = PMC
k0_all = linspace(1.0*k0*a, 1*k0*a, 40)./a;
k0_all = k0;
guessk = linspace(0.1*k0, k0*2, 200);
k_all = zeros(size(k0_all));

% pml options
h_pml = 50; % nm
PML_options = [0 h_pml 600 2];

for ii = 1:length(k0_all);
    
    fprintf('\nLoop %i of %i \n', ii, length(k0_all));
    
    % run modesolver
    % try different values of guessk
    for i_guessk = 1:length(guessk)
        
        try
            [Phi_1D, k] = complexk_mode_solver_2D_PML(n,dxz,k0_all(ii),modes,guessk(i_guessk),BC,PML_options);
        catch
            fprintf('Guess k of %s * k_0 did not work\n', guessk(i_guessk)/k0);
            continue    % continue looping until a mode is found
        end
        k_all(ii) = k(1);   % assuming the 1st solved k is the mode?
        guessk = k(1);
    
    end
%     if guessk > pi/a
% %         break
%         guessk = pi/a
%     end

end
    
% -------------------------------------------------------------------------
% PLOTTING

% reshape/normalize field
field = reshape(Phi_1D(:,1), fliplr(size(n)));  % dimensions are z, x
field = field.';    % swap dimensions
field = field./max(abs(field(:)));

% plot field for one period
figure; imagesc(imag(field)); set(gca, 'Ydir','normal'); colorbar;
title('Field, imag, one cell');
colormap(jet); 
caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]); 

% plot field across several periods
n_periods = 10;
field = repmat(field,1,n_periods);
z = linspace(0, n_periods*a, size(field,2) );
plane_wave = exp(1i*k(1).*z);
plane_wave = repmat(plane_wave,size(field,1),1);
tot_field = field.*plane_wave;

% plot imag
figure; imagesc(imag(tot_field) ); set(gca, 'Ydir','normal'); colorbar;
colormap(jet);
title('total field, imag');

% plot real
figure; imagesc(real(tot_field) ); set(gca, 'Ydir','normal'); colorbar;
colormap(jet);
title('total field, real');

% plot abs
figure; imagesc(abs(tot_field) ); set(gca, 'Ydir','normal'); colorbar;
colormap(jet);
title('total field, abs');



k

% plot vs. k
figure;
plot( real(k_all)*a, k0_all*a, 'r-o'); hold on;
plot( imag(k_all)*a, k0_all*a, 'b-o');
legend('real', 'imag');
xlabel('ka'); ylabel('k_0a');
makeFigureNice();
