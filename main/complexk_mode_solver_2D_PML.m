function [Phi_all, k_all, A, B] = complexk_mode_solver_2D_PML( N, disc, k0, num_modes, guess_k, BC, PML_options, DEBUG )
% FDFD 2D complex-k mode solver
% 
% authors: bohan zhang
%
% A reworking of Jelena's FDFD solver
% Solves for modes of 2D (transverse into plane invariant) index profile
% with 1D periodicity along propagation direction
% Allows for PEC, PMC, PML boundary conditions
% 
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
%       TEMPORARY: adding a 5th option to the PML, which is to use type 1
%       or 2
%       
%   DEBUG
%       type: boolean
%       desc: optional flag, turn on to enter debugging mode
%
% OUTPUTS:
%   Phi_all
%       type: double, ny vs. nx vs. mode #
%       desc: Output field, with dimensions ( y, x, mode # )
%   k_all
%       type: double, vector
%       desc: vector of complex wavevector eigenvalues vs. mode #

% default DEBUG to off
if nargin < 8
    DEBUG = false;
end

% total length of unwrapped vectors
[ ny, nx ]  = size(N);
n_elem      = nx*ny;

% relative permittivity
er = N.^2;

% draw in PMLs
if PML_options(1) == 1
    % currently using uniaxial pml formulation
    
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
    
    % TEMPORARY CHOOSE THE PML TYPE
    if PML_options(5) == 1
        % type 1
        
        % using slide 39 of lecture 9 slides of dr. rumpfs CEM lectures
        % setup amplitude
        ay = 1 + pml_str * ( y_indx./ny_pml ).^( pml_order );
        % setup conductivity
        sigmay = ( sin( pi*y_indx./(2*ny_pml) ).^2 );
        % combine
        eta0    = 376.73031346177;                              % ohms
        pml_y   = ( ay.*( 1 + 1i * eta0 * sigmay ) ).';
        
    elseif PML_options(5) == 2
        % type 2
        
        % using polynomial strength pml
        c       = (3e8) * (1e9);                    % nm/s
        eps0    = (8.854187817e-12) * (1e-9);       % F/nm
        omega   = k0*c;                             % rad/s
        pml_y   = 1 + ( ( 1i/( omega*eps0 ) ) * pml_str * ( y_indx./ny_pml ).^( pml_order ) ).';
        
    end
    
    % draw stretched coordinate pml
    pml_y_all                           = ones( size(N) );
    pml_y_all( 1:ny_pml, : )            = repmat( flipud(pml_y), 1, nx );
    pml_y_all( end-ny_pml+1:end, : )    = repmat( pml_y, 1, nx );
    
    % stretched coordinate operator
    Sy = spdiags( pml_y_all(:), 0, n_elem, n_elem );
    
%     % fill in pmls, working with permittivity OLD
%     er( 1:ny_pml, : )           = er( 1:ny_pml, : ).*repmat( flipud(pml_y), 1, nx );
%     er( end-ny_pml+1:end, : )   = er( end-ny_pml+1:end, : ).*repmat( pml_y, 1, nx );
    
    % DEBUG plot the pml profile
    if DEBUG
        % plot imag
        figure;
        plot( y_indx, imag(pml_y), '-o' );
        xlabel('position'); ylabel('pml profile, imag component');
        title('DEBUG imag component of pml profile');
        makeFigureNice();
        % plot real
        figure;
        plot( y_indx, real(pml_y), '-o' );
        xlabel('position'); ylabel('pml profile, real component');
        title('DEBUG real component of pml profile');
        makeFigureNice();
        % plot imag of Sy
        figure;
        plot( 1:n_elem, imag(diag(Sy)) );
        xlabel('position'); ylabel('pml profile, imag component');
        title('DEBUG imag component of diag(S_y), the stretching operator');
        makeFigureNice();
    end
    
end

% % DEBUG plot new N
% figure;
% imagesc( imag(er) );
% colorbar;
% title('DEBUG imag(\epsilon_r)');
% figure;
% imagesc( real(er) );
% colorbar;
% title('DEBUG real(\epsilon_r)');

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
if BC == 1
    % PMC
    diag0( ny:ny:end)  = 0;
end
% shift diag1
diag1( 2:end ) = diag1( 1:end-1 );
% stitch together the diags
% diag_all        = [ diagBC, diag0, diag1 ];
diag_all        = [ diag0, diag1 ];
% diag_indexs     = [ -(ny-1), 0, 1 ];
diag_indexs     = [ 0, 1 ];
% make sparse matrix
Dy_f    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% stretched coordinate pml
if PML_options(1) == 1
    Dy_f = Sy * Dy_f;
end


% generate backwards Dy
diag0               = ones( n_elem, 1 );
diagm1              = -ones( n_elem, 1 );
diagm1( ny:ny:end ) = 0;                            % no need to shift due to being in lower triangle of Dy
if BC == 1
    % PMC
    diag0( 1:ny:end ) = 0;
end
% diagBC              = zeros(n_elem, 1);             % diag 0 + (ny-1)
% diagBC( 1:ny:end )  = 0;                            % again, force PEC for now
% stitch together the diags
% diag_all        = [ diagm1, diag0, diagBC ];
diag_all        = [ diagm1, diag0 ]; 
% diag_indexs     = [ -1, 0, ny-1 ];
diag_indexs     = [ -1, 0 ];
% make sparse matrix
Dy_b    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% stretched coordinate pml
if PML_options(1) == 1
    Dy_b = Sy * Dy_b;
end

% generate Dy squared
Dy2 = Dy_b*Dy_f;

% % DEBUG plot Dy f and Dy b 
% figure;
% spy( Dy_f );
% title('DEBUG d_y forward');
% figure;
% spy( Dy_b );
% title('DEBUG d_y backward');



% generate Dx forward
diag0           = -ones( n_elem, 1 );
diagP           = ones( n_elem, 1 ); 
diagBC          = ones( n_elem, 1 );                        % BLOCH boundary conditions
diag_all        = [ diagBC, diag0, diagP ];
diag_indexes    = [ -(n_elem-ny), 0, ny ];
% make sparse matrix
Dx_f    = (1/disc)*spdiags( diag_all, diag_indexes, n_elem, n_elem );


% generate Dx backward
diag0           = ones( n_elem, 1 );
diagM           = -ones( n_elem, 1 ); 
diagBC          = -ones( n_elem, 1 );                       % BLOCH boundary conditions
diag_all        = [ diagM, diag0, diagBC ];
diag_indexes    = [ -ny, 0, (n_elem-ny) ];
% make sparse matrix
Dx_b    = (1/disc)*spdiags( diag_all, diag_indexes, n_elem, n_elem );

% % DEBUG plot Dx f and Dx b 
% figure;
% spy( Dx_f );
% title('DEBUG d_x forward');
% figure;
% spy( Dx_b );
% title('DEBUG d_x backward');


% generate Dx squared
Dx2 = Dx_b*Dx_f;


% make eigenvalue eq
n2      = spdiags( er(:), 0, n_elem, n_elem );
A       = Dx2 + Dy2 + (k0^2) * n2;
B       = 1i * ( Dx_b + Dx_f );
I       = speye( n_elem, n_elem );
Z       = sparse( n_elem, n_elem );                                 % zeros
LH      = [ A, B; Z, I ];                                           % left hand side of eigeneq
RH      = [ Z, I; I, Z ];                                           % right hand side of eigeneq

% DEBUG show these
% n2_full = full(n2);
% A_full  = full(A);
% B_full  = full(B);
% C_full  = full(C);
% Z_full  = full(Z);
% LH_full = full(LH);

% solve eigs
% Phi_out is ( Ey, Ex )(:)
[Phi_all, k_all]    = eigs(LH, RH, num_modes, guess_k);
k_all               = diag(k_all);

% unwrap the field and stuff
% reshape and sort the Phis
% when they come out raw from the modesolver, Phi_all's columns are the
% eigenvectors
% The eigenvectors are wrapped by column, then row
Phi_all     = Phi_all( 1:end/2, : );                              % first remove redundant bottom half
Phi_all     = reshape( Phi_all, ny, nx, size(Phi_all, 2) );       % hopefully this is dimensions y vs. x vs. mode#


end

