function [Phi_all, k_all, A, B] = f_bloch_complexk_mode_solver_2D_PML( N, disc, k0, num_modes, guess_k, BC, pol, PML_options, DEBUG )
% FDFD 2D complex-k mode solver
% 
% authors: bohan zhang
%
% A reworking of Jelena's FDFD solver
% Solves for modes of 2D (transverse into plane invariant) index profile
% with 1D periodicity along propagation direction
% Allows for PEC, PMC, PML boundary conditions
% 
% updated coordiantes to use x for in plane transverse, and z for in plane
% propagation
%
% INPUTS:
%   N
%       type: double, matrix
%       desc: Index distribution, dimensions x (transvers) vs z
%               (longtidinal)
%   disc
%       type: double, scalar
%       desc: discretization size (spatial units are arbitrary, as long as consistent across all inputs)
%   k0
%       type: double, scalar
%       desc: free-space wavevector (1/spatial units)
%   num_modes
%       type: integer, scalar
%       desc: number of modes to solve for
%   guess_k
%       type: double, scalar
%       desc: guess value of complex k
%   BC
%       type: integer, scalar
%       desc: y boundary condition, 0 for PEC, 1 for PMC
%   pol
%       type: string,
%       desc: polarization to solve for, 'TE' or 'TM'
%   PML_opts
%       type: 1x4 array
%       desc: PML options
%               PML_options(1): PML in y direction (yes=1 or no=0)
%               PML_options(2): length of PML layer in nm
%               PML_options(3): strength of PML in the complex plane
%               PML_options(4): PML polynomial order (1, 2, 3...)
%       TEMPORARY: adding a 5th option to the PML, which is to use type 1
%       or 2 (deprecated)
%       
%   DEBUG
%       type: boolean
%       desc: optional flag, turn on to enter debugging mode
%
% OUTPUTS:
%   Phi_all
%       type: double, nx vs. nz vs. mode #
%       desc: Output field, with dimensions ( x, z, mode # )
%   k_all
%       type: double, vector
%       desc: vector of complex wavevector eigenvalues vs. mode #

% default DEBUG to off
if nargin < 9
    DEBUG = false;
end

% total length of unwrapped vectors
[ nx, nz ]  = size(N);
n_elem      = nx*nz;

% relative permittivity
er      = N.^2;
er_yee  = repelem( er, 2, 2 );     % double grid, naive repeating

% draw in PMLs
if PML_options(1) == 1
    
    % grab params
    pml_len_nm  = PML_options(2);   % length of pml in nm
    pml_str     = PML_options(3);   % strength of pml in complex plane
    pml_order   = PML_options(4);   % pml polynomial order
    
    % setup discretizations
    nx_pml = 2 * round(pml_len_nm/disc);                                             % number of discretizations that pml spans, double sampled grid
    if abs(nx_pml - round(nx_pml)) >= 1e-5
        % discretization was not integer value
        error('Integer # of discretizations did not fit into the PML');
    end
    x_indx = 1:nx_pml;
        
    % using polynomial strength pml
%     c       = (3e8) * (1e9);                    % nm/s
%     eps0    = (8.854187817e-12) * (1e-9);       % F/nm
%     omega   = k0*c;                             % rad/s
    pml_x   = (1 + 1i * pml_str * ( x_indx./nx_pml ).^( pml_order )).';
        
    
    % draw stretched coordinate pml
    pml_x_all                           = ones( 2*size(N,1), size(N,2) );
    pml_x_all( 1:nx_pml-1, : )          = repmat( flipud( pml_x(1:end-1) ), 1, nz );
    pml_x_all( end-nx_pml+1:end, : )    = repmat( pml_x, 1, nz );
    
    % stretched coordinate operator
    pml_x_all_vec   = pml_x_all(:);
    Sx_f            = spdiags( 1./pml_x_all_vec(2:2:end), 0, n_elem, n_elem );                % half step for forward Sx
    Sx_b            = spdiags( 1./pml_x_all_vec(1:2:end-1), 0, n_elem, n_elem );              % on grid for backwards Sx
    
    % DEBUG plot the pml profile
    if DEBUG
        
        % plot imag
        figure;
        plot( x_indx, 10*log10(imag(pml_x)+1), '-o' );
        xlabel('position'); ylabel('pml profile, imag component');
        title('DEBUG imag component of pml profile + 1, in dB');
        makeFigureNice();
        
        % plot real
        figure;
        plot( x_indx, real(pml_x), '-o' );
        xlabel('position'); ylabel('pml profile, real component');
        title('DEBUG real component of pml profile');
        makeFigureNice();
        
%         % plot imag of Sy_f, dB
%         figure;
%         plot( 1:n_elem, 10*log10(imag(diag(Sy_f))+1) );
%         xlabel('position'); ylabel('pml profile, imag component');
%         title('DEBUG imag component of diag(Sy_f)+1, in dB, the forward stretching operator');
%         makeFigureNice();
        
        % plot imag of Sy_f
        figure;
        plot( 1:n_elem, imag(diag(Sx_f)) );
        xlabel('position'); ylabel('pml profile, imag component');
        title('DEBUG imag component of diag(Sy_f), the forward stretching operator');
        makeFigureNice();
        
        % plot imag of Sy_b
        figure;
        plot( 1:n_elem, imag(diag(Sx_b)) );
        xlabel('position'); ylabel('pml profile, imag component');
        title('DEBUG imag component of diag(Sy_b), the backward stretching operator');
        makeFigureNice();
        
        % plot real pml_y in 2D
        figure;
        imagesc( real(pml_x_all) );
        set( gca, 'ydir', 'normal' );
        colorbar;
        title('DEBUG real pml_y in 2D');
        makeFigureNice();
        
        % plot imag pml_y in 2D
        figure;
        imagesc( 10*log10(imag(pml_x_all)+1) );
        set( gca, 'ydir', 'normal' );
        colorbar;
        title('DEBUG imag pml_y + 1 in 2D (dB)');
        makeFigureNice();
        
    end         % end debug code
    
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

% generate forward Dx
% first generate vectors corresponding to the diagonals
% where the suffix # of the diag = # of diagonals from the middle
% Default BC is TE/PEC or TM/TMC
diag0               = -ones( n_elem, 1 );           % diag middle
% diag0( 1:nx:end )   = -sqrt(2);                            % power conserving shit
diag1               = ones( n_elem, 1 );            % diag plus 1
diag1( nx:nx:end )  = 0;                            % don't carry over into the next column
diagm1              = zeros( n_elem, 1 );
% boundary condition for TE/PMC or TM/PEC
if ( BC == 1 && strcmp( pol, 'TE' ) ) || ( BC == 0 && strcmp( pol, 'TM' ) )
    diagm1( nx-1:nx:end) = 1;     % this is using image theory
end
% shift diag1
diag1( 2:end )  = diag1( 1:end-1 );
% stitch together the diags
diag_all        = [ diagm1, diag0, diag1 ];
diag_indexs     = [ -1, 0, 1 ];
% make sparse matrix
Dx_f    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% stretched coordinate pml
if PML_options(1) == 1
    Dx_f = Sx_f * Dx_f;
end


% generate forward Dx on +1/2 grid
% Default BC is TE/PEC or TM/TMC
% the other BC is not implemented yet
diag0               = -ones( n_elem, 1 );           % diag middle
diag0( nx:nx:end )  = -2;
diag1               = ones( n_elem, 1 );            % diag plus 1
diag1( nx:nx:end )  = 0;                            % don't carry over into the next column
diagm1              = zeros( n_elem, 1 );
% % boundary condition for TE/PMC or TM/PEC
% if ( BC == 1 && strcmp( pol, 'TE' ) ) || ( BC == 0 && strcmp( pol, 'TM' ) )
%     diagm1( nx-1:nx:end) = 1;     % this is using image theory
% end
% shift diag1
diag1( 2:end )  = diag1( 1:end-1 );
% stitch together the diags
diag_all        = [ diagm1, diag0, diag1 ];
diag_indexs     = [ -1, 0, 1 ];
% make sparse matrix
Dx_f_plushalf    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% stretched coordinate pml
if PML_options(1) == 1
    Dx_f_plushalf = Sx_b * Dx_f_plushalf;
end


% % Forward Dx, with opposite boundary condition
% diag0               = -ones( n_elem, 1 );           % diag middle
% diag1               = ones( n_elem, 1 );            % diag plus 1
% diag1( ny:ny:end )  = 0;                            % field outside domain = 0?
% diagm1              = zeros( n_elem, 1 );
% diagm1( ny-1:ny:end) = 1;     % this is using image theory
% % boundary condition for TE/PMC or TM/PEC
% if ( BC == 1 && strcmp( pol, 'TE' ) ) || ( BC == 0 && strcmp( pol, 'TM' ) )
%     % is this right?
%     diagm1              = zeros( n_elem, 1 );
% end
% % shift diag1
% diag1( 2:end )  = diag1( 1:end-1 );
% % stitch together the diags
% diag_all        = [ diagm1, diag0, diag1 ];
% diag_indexs     = [ -1, 0, 1 ];
% % make sparse matrix
% Dy_f_opposite_bc    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% % stretched coordinate pml
% if PML_options(1) == 1
%     Dy_f_opposite_bc = Sy_f * Dy_f_opposite_bc;
% end


% generate backwards Dx
diag0               = ones( n_elem, 1 );
% diag0( nx:nx:end )   = sqrt(2);                            % power conserving shit
diagm1              = -ones( n_elem, 1 );
diagm1( nx:nx:end ) = 0;                            % no need to shift due to being in lower triangle of Dx
diag1               = zeros( n_elem, 1 );
% boundary condition for TE/PMC or TM/PEC
if ( BC == 1 && strcmp( pol, 'TE' ) ) || ( BC == 0 && strcmp( pol, 'TM' ) )
    diag1( 1:nx:end ) = -1;    % this is using image theory
end
% shift diag1
diag1( 2:end )  = diag1( 1:end-1 ); 
% stitch together the diags
diag_all        = [ diagm1, diag0, diag1 ]; 
diag_indexs     = [ -1, 0, 1 ];
% make sparse matrix
Dx_b    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% stretched coordinate pml
if PML_options(1) == 1
    Dx_b = Sx_b * Dx_b;
end

% % generate backwards Dy, with opposite BC
% diag0               = ones( n_elem, 1 );
% diagm1              = -ones( n_elem, 1 );
% diagm1( ny:ny:end ) = 0;                            % no need to shift due to being in lower triangle of Dy
% diag1               = zeros( n_elem, 1 );
% diag1( 1:ny:end ) = -1;    % this is using image theory
% % boundary condition for TE/PMC or TM/PEC
% if ( BC == 1 && strcmp( pol, 'TE' ) ) || ( BC == 0 && strcmp( pol, 'TM' ) )
% %     diag0( ny:ny:end)  = 0;
%     diag1               = zeros( n_elem, 1 );
% end
% % shift diag1
% diag1( 2:end )  = diag1( 1:end-1 ); 
% % stitch together the diags
% diag_all        = [ diagm1, diag0, diag1 ]; 
% diag_indexs     = [ -1, 0, 1 ];
% % make sparse matrix
% Dy_b_opposite_bc    = (1/disc)*spdiags( diag_all, diag_indexs, n_elem, n_elem );
% % stretched coordinate pml
% if PML_options(1) == 1
%     Dy_b_opposite_bc = Sy_b * Dy_b_opposite_bc;
% end

% generate Dx squared
Dx2 = Dx_b*Dx_f;

% generate Dx center
Dx_center = (Dx_f + Dx_b)/2;

% % DEBUG plot Dy f and Dy b 
% figure;
% spy( Dy_f );
% title('DEBUG d_y forward');
% figure;
% spy( Dy_b );
% title('DEBUG d_y backward');

% generate Dz forward
diag0           = -ones( n_elem, 1 );
diagP           = ones( n_elem, 1 ); 
diagBC          = ones( n_elem, 1 );                        % BLOCH boundary conditions
diag_all        = [ diagBC, diag0, diagP ];
diag_indexes    = [ -(n_elem-nx), 0, nx ];
% make sparse matrix
Dz_f    = (1/disc)*spdiags( diag_all, diag_indexes, n_elem, n_elem );

% generate Dx backward
diag0           = ones( n_elem, 1 );
diagM           = -ones( n_elem, 1 ); 
diagBC          = -ones( n_elem, 1 );                       % BLOCH boundary conditions
diag_all        = [ diagM, diag0, diagBC ];
diag_indexes    = [ -nx, 0, (n_elem-nx) ];
% make sparse matrix
Dz_b    = (1/disc)*spdiags( diag_all, diag_indexes, n_elem, n_elem );

% % DEBUG plot Dx f and Dx b 
% figure;
% spy( Dx_f );
% title('DEBUG d_x forward');
% figure;
% spy( Dx_b );
% title('DEBUG d_x backward');

% generate Dz squared
Dz2 = Dz_b*Dz_f;

% generate Dz center
Dz_center = (Dz_f + Dz_b)/2;

% take derivatives of permittivity
% wait what about pmls
er_extended_z   = [ er(:,end), er, er(:,1) ];                               % CAREFUL THIS HAS TO BE PERIODIC
dz_er           = ( er_extended_z( :, 3:end ) - er_extended_z( :, 1:end-2 ) )./(2*disc);
er_extended_x   = [ er(1,:); er; er(end,:) ];
dx_er           = ( er_extended_x( 3:end, : ) - er_extended_x( 1:end-2, : ) )./(2*disc);

dz_er_center = Dz_center * er(:);
dx_er_center = Dx_center * er(:);

% % DEBUG plotting the two derivatives of epsilon
% figure;
% plot( 1:length(dx_er(:)), dx_er(:) ); hold on;
% plot( 1:length(dx_er_center), dx_er_center );
% makeFigureNice();

% epsilon = n^2 operator
n2      = spdiags( er(:), 0, n_elem, n_elem );
n2_inv  = spdiags( 1./er(:),  0, n_elem, n_elem );

% identity and zero operators
I       = speye( n_elem, n_elem );
Z       = sparse( n_elem, n_elem );

% make eigenvalue eq



if strcmp( pol, 'TE' )
    
    % TE wave equation
    A = Dx2 + Dz2 + (k0^2) * n2;
    B = 1i * ( Dz_b + Dz_f );
    C = -1*speye( n_elem, n_elem );
    
elseif strcmp( pol, 'TM' )
    
    er_inv  = 1./er;
    n2_inv  = spdiags( er_inv(:), 0, n_elem, n_elem );
    
    % TM wave equation
%     Dz      = ( Dz_b + Dz_f )./2;
%     Dx      = ( Dx_b + Dx_f )./2;
%     Dy_opp  = ( Dy_b_opposite_bc + Dy_f_opposite_bc )./2;
%     Dx_n2   = spdiags( Dx * er(:), 0, n_elem, n_elem );
%     Dx_n2_inv   = spdiags( Dx * er_inv(:), 0, n_elem, n_elem );
%     Dy_n2_old   = spdiags( Dy * er(:), 0, n_elem, n_elem );
%     Dx_n2   = spdiags( dx_er(:), 0, n_elem, n_elem );
%     Dy_n2   = spdiags( dx_er(:), 0, n_elem, n_elem );
%     A       = Dx2 + Dy2 + (k0^2)*n2 - n2_inv*Dx_n2*Dx - n2_inv*Dy_n2*Dy;
%     B       = 2*1i*Dx - 1i*n2_inv*Dx_n2;
% 
%     % alternative formulation that doesn't require expansion across y
%     A       = Dx2 + (k0^2)*n2 - n2_inv*Dx_n2*Dx + n2*Dy*n2_inv*Dy;
%     B       = 2*1i*Dx - 1i*n2_inv*Dx_n2;
    
%     % alternative formulation that uses dX(1/epsilon)
%     Dx_f_n2_inv = 
%     A = n2*Dx_n2_inv*Dx + Dx2 + n2*Dy*n2_inv*Dy + (k0^2)*n2;
%     B = n2*Dx_n2_inv*1i + 2i*Dx;
    
%     % alternative formulation that uses dX(1/epsilon) and opposite bc
%     A = n2*Dx_n2_inv*Dx + Dx2 + n2*Dy_opp*n2_inv*Dy + (k0^2)*n2;
%     B = n2*Dx_n2_inv*1i + 2i*Dx;
%     
%     % alternative formulation that uses dX(1/epsilon) and expands dy
%     A = n2*Dx_n2_inv*Dx + Dx2 + Dy2 - n2_inv*Dy_n2*Dy + (k0^2)*n2;
%     B = n2*Dx_n2_inv*1i + 2i*Dx;

%     % alternative formulation that is some bullshit
%     A = n2*Dx*n2_inv*Dx + Dx2 + n2*Dy*n2_inv*Dy + (k0^2)*n2;
%     B = n2*Dx*n2_inv*1i + 2i*Dx;

%     % alternative formulation that uses dX(1/epsilon) and Dyf and Dyb
%     % this one was working the best.
%     Dz_n2_inv   = spdiags( Dz * er_inv(:), 0, n_elem, n_elem );
%     A = n2*Dz_n2_inv*Dz + Dz2 + n2*Dx_b*n2_inv*Dx_f + (k0^2)*n2;
%     B = n2*Dz_n2_inv*1i + 2i*Dz;
%     C = -1*speye( n_elem, n_elem );
    
    % latest latest version
%     er_m12_p1           = [ er_yee( 2:2:end, 3:2:end ), er_yee( 2:2:end, 1 ) ];     % er(m+1/2, p+1)
%     er_m12_p            = er_yee( 2:2:end, 1:2:end );                               % er(m+1/2, p)
%     er_m_p12            = er_yee( 1:2:end, 2:2:end );                               % er(m, p+1/2)
%     er_m12_p1_inv       = 1./er_m12_p1;
%     er_m12_p_inv        = 1./er_m12_p;
%     er_m_p12_inv        = 1./er_m_p12;
%     dz_f_er_m12_p_inv   = Dz_f * er_m12_p_inv(:);
%     
%     % in matrix format
%     er_m12_p1_inv_sparse        = spdiags( er_m12_p1_inv(:), 0, n_elem, n_elem );
%     er_m12_p_inv_sparse         = spdiags( er_m12_p_inv(:), 0, n_elem, n_elem );
%     dz_f_er_m12_p_inv_sparse    = spdiags( dz_f_er_m12_p_inv(:), 0, n_elem, n_elem );
%     er_m_p12_inv_sparse         = spdiags( er_m_p12_inv(:), 0, n_elem, n_elem );
% 
%     % build operators
% %     A = er_m12_p1_inv_sparse*Dz_f*Dz_b + dz_f_er_m12_p_inv_sparse*Dz_b + Dx_f*er_m_p12_inv_sparse*Dx_b + k0.^2;
%     A = er_m12_p_inv_sparse*Dz_f*Dz_b + dz_f_er_m12_p_inv_sparse*Dz_b + Dx_f*er_m_p12_inv_sparse*Dx_b + k0.^2;
%     B = -1i*er_m12_p_inv_sparse*Dz_b - 1i*er_m12_p1_inv_sparse*Dz_f - 1i*dz_f_er_m12_p_inv_sparse;
%     B = -1i*er_m12_p_inv_sparse*Dz_b - 1i*er_m12_p_inv_sparse*Dz_f - 1i*dz_f_er_m12_p_inv_sparse;
%     C = -er_m12_p_inv_sparse;

    A = n2 * Dz_f * n2_inv * Dz_b + n2 * Dx_f_plushalf * n2_inv * Dx_b + (k0.^2) * n2;
    B = -1i * ( Dz_b + n2 * Dz_f * n2_inv );
    C = -I;

else
    error('Polarization input must be ''TE'' or ''TM'' ');
end                              
LH      = [ A, B; Z, I ];                                           % left hand side of eigeneq
RH      = [ Z, -C; I, Z ];                                          % right hand side of eigeneq

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
Phi_all     = reshape( Phi_all, nx, nz, size(Phi_all, 2) );       % hopefully this is dimensions x vs. z vs. mode#


end

