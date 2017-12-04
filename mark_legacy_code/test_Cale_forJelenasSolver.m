% editing by bz
% all units nm

clear; close all; clc;

% dimensions?
a = 0.330;
% a = 0.5;
b = 0.600;
c = 0.200;

% indices
n1 = 1.45;
n2 = 3.5;

% not sure what this is
zf = a;
xf = 3*b;

% discretization
dx = 0.01; dz = 0.01;

% index refraction
n = n1*ones(round(xf/dx),round(zf/dz));     % dimensions (x, z)
n(round(b/dx)+1:round(2*b/dx),:) = n2;
n(round((1.5*b-.5*c)/dx):round((1.5*b+.5*c)/dx),1:round((c)/dz)) = n1; 

% display
figure; imagesc(n); set(gca,'Ydir','normal'); colorbar;
xlabel('z'); ylabel('x');

% inputs
guessk = 0.8*pi/(a*1e3);  %switch to nm, why lol
k0 = 2*pi/1550;
d = dx*1e3;
modes = 10;
BC = 0;

% pml options
h_pml = 100; 
PML_options = [1 h_pml 600 2];

% run modesolver
[Phi_1D, k] = complexk_mode_solver_2D_PML(n,d,k0,modes,guessk,BC,PML_options);

% plot field
field = reshape(Phi_1D(:,1), fliplr(size(n)));

figure; imagesc(imag(field.')); set(gca, 'Ydir','normal'); colorbar;
% colormap('redbluehilight'); 
caxis([-max(max(abs(imag(field)))) max(max(abs(imag(field))))]); 

k