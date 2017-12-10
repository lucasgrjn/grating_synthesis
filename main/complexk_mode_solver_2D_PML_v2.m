function [Phi_1D, k] = complexk_mode_solver_2D_PML_v2( N, disc, k0, num_modes, guess_k, BC, PML_options )
% FDFD 2D complex-k mode solver
% 
% authors: bohan zhang
%
% A reworking of Jelena's FDFD solver

% Jelena's FDFD 2D complex-k mode solver
% Version 3
% May 15, 2014
% Includes: Perfect Electric or Magnetic Boundary Conditions, Periodic
% Boundary Conditions, Pefectly Matched Layers (PMLs)
%
% INPUTS:
%   N
%       type: double, matrix
%       desc: Index distribution
%   disc
%       type: double, scalar
%       desc: discretization size (nm)
%   k0
%       type: double, scalar
%       desc: free-space wavevector (1/nm)
%   num_modes
%       type: integer, scalar
%       desc: number of modes to solve for
%   guess_k
%       type: double, scalar
%       desc: guess value of complex k
%   BC
%       type: integer, scalar
%       desc: y boundary condition, 0 for PEC, 1 for PMC
%   PML_opts
%       type: 1x4 array
%       desc: PML options
%               PML_options(1): PML in y direction (yes=1 or no=0)
%               PML_options(2): length of PML layer in nm
%               PML_options(3): strength of PML in the complex plane
%               PML_options(4): PML polynomial order (1, 2, 3...)
%
% OUTPUTS:
% Phi_1D: matrix of all Phi-field solutions (E = e^(jkx)*Phi) in 1D form for desired structure
% k: vector of complex wavevector values calulated

% total length of unwrapped vectors
[ ny, nx ]  = size(N);
n_elem      = nx*ny;


% draw in PMLs
if PML_options(1) == 1
   
    % grab params
    pml_len_nm  = PML_options(2);   % length of pml in nm
    pml_str     = PML_options(3);   % strength of pml in complex plane
    pml_order   = PML_options(4);   % pml polynomial order
    
    % setup discretizations
    ny_pml = pml_len_nm/disc;                      % number of discretizations that pml spans
    if abs(ny_pml - round(ny_pml)) >= 1e-5
        % discretization was not integer value
        error('Integer # of discretizations did not fit into the PML');
    end
    y_indx = 1:ny_pml;
    
    % using slide 39 of lecture 9 slides of dr. rumpfs CEM lectures
    % setup amplitude
    ay = 1 + pml_str * ( y_indx./ny_pml ).^( pml_order );
    % setup conductivity
    sigmay = ( sin( pi*y_indx./(2*ny_pml) ).^2 );
    % combine
    eta0    = 376.73031346177;                              % ohms
    pml_y   = ( ay.*( 1 + 1i * eta0 * sigmay ) ).';
    
    % fill in pmls
    N( 1:ny_pml, : )            = N( 1:ny_pml, : ).*repmat( flipud(pml_y), 1, nx );
    N( end-ny_pml+1:end, : )    = N( end-ny_pml+1:end, : ).*repmat( pml_y, 1, nx );
    
end

% DEBUG plot new N
figure;
imagesc( imag(N) );
colorbar;
title('DEBUG imag(N)');
figure;
imagesc( real(N) );
colorbar;
title('DEBUG real(N)');

% the rule of spdiags is:
% In this syntax, if a column of B is longer than the diagonal it is
% replacing, and m >= n, (m = num of rows, n = num of cols in the SPARSE matrix)
% spdiags takes elements of super-diagonals from the lower part of the column of B, and 
% elements of sub-diagonals from the upper part of the column of B. 
% However, if m < n , then super-diagonals are from the upper part of the column of B, 
% and sub-diagonals from the lower part.
% since m = n for us, always, then we are always replacing with the lower
% part of the column

% generate forward Dy
% first generate vectors corresponding to the diagonals
% where the suffix # of the diag = # of diagonals from the middle
diag0               = -ones( n_elem, 1 );           % diag middle
diag1               = ones( n_elem, 1 );            % diag plus 1
diag1( ny:ny:end )  = 0;

% TEMP: this boundary condition will later be either PMC (val 1) or PEC
% (val 0 )
diagBC              = zeros(n_elem, 1);             % diag 0 - (ny-1)
diagBC( 1:ny:end )  = 0;                            % force PEC for now

% shift diag1
diag1( 2:end ) = diag1( 1:end-1 );

% stitch together the diags
diag_all        = [ diagBC, diag0, diag1 ];
diag_indexs     = [ -(ny-1), 0, 1 ];

% make sparse matrix
Dy_f    = spdiags( diag_all, diag_indexs, n_elem, n_elem );

% DEBUG show Dy_f
% full(Dy_f)


% generate backwards Dy
diag0               = ones( n_elem, 1 );
diagm1              = -ones( n_elem, 1 );
diagm1( ny:ny:end ) = 0;                            % no need to shift due to being in lower triangle of Dy
diagBC              = zeros(n_elem, 1);             % diag 0 + (ny-1)
diagBC( 1:ny:end )  = 0;                            % again, force PEC for now

% stitch together the diags
diag_all        = [ diagm1, diag0, diagBC ];
diag_indexs     = [ -1, 0, ny-1 ];

% make sparse matrix
Dy_b    = spdiags( diag_all, diag_indexs, n_elem, n_elem );

% % DEBUG show Dy_b
% full(Dy_b)

% multiply Dy_f and Dy_b
Dy2 = Dy_b*Dy_f;

% % DEBUG show Dy2
% full(Dy2)


% generate Dx forward
diag0   = -ones( n_elem, 1 );
diagP   = ones( n_elem, 1 ); 
diagBC  = ones( n_elem, 1 );                        % BLOCH boundary conditions
diag_all        = [ diagBC, diag0, diagP ];
diag_indexes    = [ -(n_elem-ny+1), 0, ny-1 ];

% make sparse matrix
Dx_f    = spdiags( diag_all, diag_indexes, n_elem, n_elem );

% % DEBUG show Dx_f
% full(Dx_f)


% generate Dx backward
diag0   = ones( n_elem, 1 );
diagM   = -ones( n_elem, 1 ); 
diagBC  = -ones( n_elem, 1 );                       % BLOCH boundary conditions
diag_all        = [ diagM, diag0, diagBC ];
diag_indexes    = [ -(ny-1), 0, (n_elem-ny+1) ];

% make sparse matrix
Dx_b    = spdiags( diag_all, diag_indexes, n_elem, n_elem );

% % DEBUG show Dx_b
% full(Dx_b)


% generate Dx squared
Dx2 = Dx_b*Dx_f;

% % DEBUG show Dx2
% full(Dx2)

% make eigenvalue eq
n2      = spdiags( N(:).^2, 0, n_elem, n_elem );
A       = Dx2 + Dy2 + (k0^2) * n2;
B       = 1i * ( Dx_b + Dx_f );
C       = -speye( n_elem, n_elem );
Z       = speye( n_elem, n_elem );                                  % zeros
LH      = [ A, B; Z, -C ];                                          % left hand side of eigeneq
RH      = [ Z, -C; -C, Z ];                                         % right hand side of eigeneq

% DEBUG show these
n2_full = full(n2);
A_full  = full(A);
B_full  = full(B);
C_full  = full(C);
Z_full  = full(Z);
LH_full = full(LH);

% solve eigs
[Phi_out, k_out]    = eigs(LH, RH, num_modes, guess_k);
k_unsorted          = diag(k_out);

% TEMP dummy code
Phi_1D = [];
k = [];

% -------------------------------------------------------------------------
% OLD CODE

% % Check function input
% if nargin ~= 7
%     error('ERROR: Incorrect number of function inputs.');
% end
% 
% % Perfect magnetic and perfectly matched boundaries
% PMC = BC;
% PML = PML_options(1);
% if PML==1
%     PMC = 0;
%     if size(PML_options) ~= 4
%         error('ERROR: PML turned on but incorrect number of PML option inputs.');
%     end
%     PML_length = PML_options(2);
%     PML_strength = PML_options(3);
%     PML_order = PML_options(4);
% end
% 
% % Compute necessary variables and setup structure vectors
% number_x = size(n,2);
% number_y = size(n,1);
% number_xy = number_x*number_y;
% width_x = number_x*d;
% width_y = number_y*d;
% a = width_x;
% x = zeros(number_x,1);
% x(1) = 0;
% for j = 2:number_x;
%     x(j) = x(j-1) + d;
% end
% y = zeros(number_y,1);
% y(1) = 0;
% for j = 2:number_y;
%     y(j) = y(j-1) + d;
% end
% if PML == 1
%     PML_points = round(PML_length/d);
%     PML_mult = PML_strength/(PML_points^PML_order);
%     j = PML_points;
%     for k=1:PML_points,
%         y(k) = y(k) - PML_mult*(j^PML_order)*1i;
%         j = j-1;
%     end
%     j = 1;
%     for k=(number_y-PML_points+1):number_y,
%         y(k) = y(k) + PML_mult*(j^PML_order)*1i;
%         j = j+1;
%     end
% end
% % Create y half space vector
% yhalf = zeros(1,number_y+1);
% yhalf(2:number_y) = y(1:number_y-1)+(y(2:number_y)-y(1:number_y-1))/2;
% yhalf(1) = yhalf(2)+(yhalf(2)-yhalf(3));
% yhalf(number_y+1) = yhalf(number_y)+(yhalf(number_y)-yhalf(number_y-1));
% 
% % % DEBUGGING
% % fprintf('Time to set up A matrix:\n');
% % tic;
% 
% % Setup A matrix
% A = sparse(number_xy);
% for k = 1:number_xy
%     k_y = floor((k+number_x-1)/number_x);
%     k_x = mod(k+number_x,number_x);
%     if k_x == 0
%         k_x = number_x;
%     end
%     if (k_y==1) && (k_x==1)
%         % first row first col
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k+number_x) = 1/d^2;
%         % periodic boundary condition
%         A(k,k+number_x-1) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k+number_x) = 2/d^2;
%         end
%     elseif (k_y==number_y) && (k_x==number_x)
%         % last row last col
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k-number_x) = 1/d^2;
%         % periodic boundary condition
%         A(k,k-number_x+1) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k-number_x) = 2/d^2;
%         end
%     elseif (k_y==number_y) && (k_x==1)
%         % last row first col
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k-number_x) = 1/d^2;
%         % periodic boundary condition
%         A(k,k+number_x-1) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k-number_x) = 2/d^2;
%         end
%     elseif (k_y==1) && (k_x==number_x)
%         % first row last col
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k+number_x) = 1/d^2;
%         % periodic boundary condition
%         A(k,k-number_x+1) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k+number_x) = 2/d^2;
%         end
%     elseif k_x==1
%         % first col
%         A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
%         A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
%         % periodic boundary condition first col
%         A(k,k+number_x-1) = 1/d^2;
%     elseif k_x==number_x
%         % last col
%         A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
%         A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
%         % periodic boundary condition last col
%         A(k,k-number_x+1) = 1/d^2;
%     elseif k_y==1
%         % first row
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k+number_x) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k+number_x) = 2/d^2;
%         end
%     elseif k_y==number_y
%         % last row
%         A(k,k) = -4/d^2 + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k-number_x) = 1/d^2;
%         % perfect magnetic boundary
%         if (PMC == 1)
%             A(k,k-number_x) = 2/d^2;
%         end
%     else
%         A(k,k) = -2/d^2 - 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y))) - 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y))) + (k0*n(k_y,k_x))^2;
%         A(k,k+1) = 1/d^2;
%         A(k,k-1) = 1/d^2;
%         A(k,k+number_x) = 1/((y(k_y+1)-y(k_y))*(yhalf(k_y+1)-yhalf(k_y)));
%         A(k,k-number_x) = 1/((y(k_y)-y(k_y-1))*(yhalf(k_y+1)-yhalf(k_y)));
%     end
% end
% 
% % toc;
% % % DEBUGGING
% % fprintf('Time to set up B matrix:\n');
% % tic;
% 
% % Setup B matrix
% B = sparse(number_xy);
% for k = 1:number_xy
%     k_y = floor((k+number_x-1)/number_x);
%     k_x = mod(k+number_x,number_x);
%     if k_x == 0
%         k_x = number_x;
%     end
%     if (k_y==1) && (k_x==1)
%         % first row first col
%         B(k,k+1) = 1i*1/d;
%         % periodic boundary condition first col
%         B(k,k+number_x-1) = -1i*1/d;
%     elseif (k_y==number_y) && (k_x==number_x)
%         % last row last col
%         B(k,k-1) = -1i*1/d;
%         % periodic boundary condition last col
%         B(k,k-number_x+1) = 1i*1/d;
%     elseif (k_y==number_y) && (k_x==1)
%         % last row first col
%         B(k,k+1) = 1i*1/d;
%         % periodic boundary condition first col
%         B(k,k+number_x-1) = -1i*1/d;
%     elseif (k_y==1) && (k_x==number_x)
%         % first row last col
%         B(k,k-1) = -1i*1/d;
%         % periodic boundary condition last col
%         B(k,k-number_x+1) = 1i*1/d;
%     elseif k_x==1
%         % first col
%         B(k,k+1) = 1i*1/d;
%         % periodic boundary condition first col
%         B(k,k+number_x-1) = -1i*1/d;
%     elseif k_x==number_x
%         % last col
%         B(k,k-1) = -1i*1/d;
%         % periodic boundary condition last col
%         B(k,k-number_x+1) = 1i*1/d;
%     elseif k_y==1
%         % first row
%         B(k,k+1) = 1i*1/d;
%         B(k,k-1) = -1i*1/d;
%     elseif k_y==number_y
%         % last row
%         B(k,k+1) = 1i*1/d;
%         B(k,k-1) = -1i*1/d;
%     else
%         B(k,k+1) = 1i*1/d;
%         B(k,k-1) = -1i*1/d;
%     end
% end
% 
% % toc;
% % % DEBUGGING
% % fprintf('Time to solve eigs:\n');
% % tic;
% 
% % Setup C matrix
% C = -1.*speye(number_x*number_y);
% 
% % Linearize Eigen problem and solve
% L1 = [A B; sparse(number_x*number_y,number_x*number_y) speye(number_x*number_y)];
% L2 = [sparse(number_x*number_y,number_x*number_y) -C; speye(number_x*number_y) sparse(number_x*number_y,number_x*number_y)];
% [Phi_out, k_out] = eigs(L1, L2, modes, guessk);
% k_unsorted = diag(k_out);
% 
% % toc;
% 
% % Sort outputs for 0<k<pi/a
% jjj = 1;
% for jj = 1:length(k_unsorted)
%     if (real(k_unsorted(jj)) <= 1.03*pi/a) && (real(k_unsorted(jj))>=0)
%         k(jjj,1) = k_unsorted(jj);
%         Phi_sorted(:,jjj) = Phi_out(:,jj);
%         jjj = jjj+1;
%     end
% end
% Phi_1D = Phi_sorted(1:length(Phi_sorted)/2,:);

end

