function [Phi_1D, k] = complexk_mode_solver_2D_PML(n,d,k0,modes,guessk,BC,PML_options)
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


% Check function input
if nargin ~= 7
    error('ERROR: Incorrect number of function inputs.');
end

% Perfect magnetic and perfectly matched boundaries
PMC = BC;
PML = PML_options(1);
if PML==1
    PMC = 0;
    if size(PML_options) ~= 4
        error('ERROR: PML turned on but incorrect number of PML option inputs.');
    end
    PML_length = PML_options(2);
    PML_strength = PML_options(3);
    PML_order = PML_options(4);
end

% Compute necessary variables and setup structure vectors
number_x = size(n,2);
number_y = size(n,1);
number_xy = number_x*number_y;
width_x = number_x*d;
width_y = number_y*d;
a = width_x;
x = zeros(number_x,1);
x(1) = 0;
for j = 2:number_x;
    x(j) = x(j-1) + d;
end
y = zeros(number_y,1);
y(1) = 0;
for j = 2:number_y;
    y(j) = y(j-1) + d;
end
if PML == 1
    PML_points = round(PML_length/d);
    PML_mult = PML_strength/(PML_points^PML_order);
    j = PML_points;
    for k=1:PML_points,
        y(k) = y(k) - PML_mult*(j^PML_order)*1i;
        j = j-1;
    end
    j = 1;
    for k=(number_y-PML_points+1):number_y,
        y(k) = y(k) + PML_mult*(j^PML_order)*1i;
        j = j+1;
    end
end
% Create y half space vector
yhalf = zeros(1,number_y+1);
yhalf(2:number_y) = y(1:number_y-1)+(y(2:number_y)-y(1:number_y-1))/2;
yhalf(1) = yhalf(2)+(yhalf(2)-yhalf(3));
yhalf(number_y+1) = yhalf(number_y)+(yhalf(number_y)-yhalf(number_y-1));

% % DEBUGGING
% fprintf('Time to set up A matrix:\n');
% tic;

% Setup A matrix
A = sparse(number_xy);
for k = 1:number_xy
    k_y = floor((k+number_x-1)/number_x);
    k_x = mod(k+number_x,number_x);
    if k_x == 0
        k_x = number_x;
    end
    if (k_y==1) && (k_x==1)
        % first row first col
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k+number_x) = 1/d^2;
        % periodic boundary condition
        A(k,k+number_x-1) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k+number_x) = 2/d^2;
        end
    elseif (k_y==number_y) && (k_x==number_x)
        % last row last col
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k-1) = 1/d^2;
        A(k,k-number_x) = 1/d^2;
        % periodic boundary condition
        A(k,k-number_x+1) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k-number_x) = 2/d^2;
        end
    elseif (k_y==number_y) && (k_x==1)
        % last row first col
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k-number_x) = 1/d^2;
        % periodic boundary condition
        A(k,k+number_x-1) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k-number_x) = 2/d^2;
        end
    elseif (k_y==1) && (k_x==number_x)
        % first row last col
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k-1) = 1/d^2;
        A(k,k+number_x) = 1/d^2;
        % periodic boundary condition
        A(k,k-number_x+1) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k+number_x) = 2/d^2;
        end
    elseif k_x==1
        % first col
        A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
        A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
        % periodic boundary condition first col
        A(k,k+number_x-1) = 1/d^2;
    elseif k_x==number_x
        % last col
        A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
        A(k,k-1) = 1/d^2;
        A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
        A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
        % periodic boundary condition last col
        A(k,k-number_x+1) = 1/d^2;
    elseif k_y==1
        % first row
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k-1) = 1/d^2;
        A(k,k+number_x) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k+number_x) = 2/d^2;
        end
    elseif k_y==number_y
        % last row
        A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k-1) = 1/d^2;
        A(k,k-number_x) = 1/d^2;
        % perfect magnetic boundary
        if (PMC == 1)
            A(k,k-number_x) = 2/d^2;
        end
    else
        A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
        A(k,k+1) = 1/d^2;
        A(k,k-1) = 1/d^2;
        A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
        A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
    end
end

% toc;
% % DEBUGGING
% fprintf('Time to set up B matrix:\n');
% tic;

% Setup B matrix
B = sparse(number_xy);
for k = 1:number_xy
    k_y = floor((k+number_x-1)/number_x);
    k_x = mod(k+number_x,number_x);
    if k_x == 0
        k_x = number_x;
    end
    if (k_y==1) && (k_x==1)
        % first row first col
        B(k,k+1) = 1i*1/d;
        % periodic boundary condition first col
        B(k,k+number_x-1) = -1i*1/d;
    elseif (k_y==number_y) && (k_x==number_x)
        % last row last col
        B(k,k-1) = -1i*1/d;
        % periodic boundary condition last col
        B(k,k-number_x+1) = 1i*1/d;
    elseif (k_y==number_y) && (k_x==1)
        % last row first col
        B(k,k+1) = 1i*1/d;
        % periodic boundary condition first col
        B(k,k+number_x-1) = -1i*1/d;
    elseif (k_y==1) && (k_x==number_x)
        % first row last col
        B(k,k-1) = -1i*1/d;
        % periodic boundary condition last col
        B(k,k-number_x+1) = 1i*1/d;
    elseif k_x==1
        % first col
        B(k,k+1) = 1i*1/d;
        % periodic boundary condition first col
        B(k,k+number_x-1) = -1i*1/d;
    elseif k_x==number_x
        % last col
        B(k,k-1) = -1i*1/d;
        % periodic boundary condition last col
        B(k,k-number_x+1) = 1i*1/d;
    elseif k_y==1
        % first row
        B(k,k+1) = 1i*1/d;
        B(k,k-1) = -1i*1/d;
    elseif k_y==number_y
        % last row
        B(k,k+1) = 1i*1/d;
        B(k,k-1) = -1i*1/d;
    else
        B(k,k+1) = 1i*1/d;
        B(k,k-1) = -1i*1/d;
    end
end

% toc;
% % DEBUGGING
% fprintf('Time to solve eigs:\n');
% tic;

% Setup C matrix
C = -1.*speye(number_x*number_y);

% Linearize Eigen problem and solve
L1 = [A B; sparse(number_x*number_y,number_x*number_y) speye(number_x*number_y)];
L2 = [sparse(number_x*number_y,number_x*number_y) -C; speye(number_x*number_y) sparse(number_x*number_y,number_x*number_y)];
[Phi_out, k_out] = eigs(L1, L2, modes, guessk);
k_unsorted = diag(k_out);

% toc;

% Sort outputs for 0<k<pi/a
jjj = 1;
for jj = 1:length(k_unsorted)
    if (real(k_unsorted(jj)) <= 1.03*pi/a) && (real(k_unsorted(jj))>=0)
        k(jjj,1) = k_unsorted(jj);
        Phi_sorted(:,jjj) = Phi_out(:,jj);
        jjj = jjj+1;
    end
end
Phi_1D = Phi_sorted(1:length(Phi_sorted)/2,:);
end

