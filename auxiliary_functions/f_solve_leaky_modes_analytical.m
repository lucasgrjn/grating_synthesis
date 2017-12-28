function [ beta ] = f_solve_leaky_modes_analytical( a, b, n0, n1, lambda0 )
% CURRENTLY DEPRECATED... since matlab cannot minimize a complex valued
% function (where the x axis is complex beta)
%
% Solves for TE analytical solutions to symmetric leaky mode waveguide
%
% Using formulation from paper:
% j hu c menyuk - understanding leaky modes slab waveguide revisited
%
% Leaky mode waveguide index profile is defined as:
%   
%   n0, upper cladding towards infinity
%   -------------------
%   n1, width b
%   ------------------- 
%   n0, width a
%   -------------------
%   n1, width b
%   -------------------
%   n0, bottom cladding towards neg. infinity
%
%   in plane transverse direction = x
%   direction of propagation = z
%
%
% authors: bohan zhang
%
% Inputs:
%   a
%       type: double, scalar
%       desc: width of waveguide core with index n0, units of nm
%   b
%       type: double, scalar
%       desc: width of waveguide inter-cladding layer with index n1, units
%             of nm
%   n0
%       type: double, scalar
%       desc: index of waveguide core/outer cladding. Should be > n1
%   n1
%       type: double, scalar
%       desc: index of waveguide inter-cladding layer. Should be < n0
%   lambda0
%       type: double, scalar
%       desc: free space wavelength in nm

% define constants
k0 = 2*pi/lambda0;  % units rad/nm

% For first test, plot the dispersion relation

% pick range of beta's to solve for
neff_range  = linspace( n1*1.0001, n0*0.9999, 1000 );
beta_range  = (2*pi/lambda0).*neff_range;
% kx          = sqrt( (k0^2)*(n0^2) - beta_range.^2 );
% alpha       = sqrt( beta_range.^2 - (k0^2)*(n1^2) );

% dispersion equation
y_lhs   = lhs( beta_range, a, n0, k0 );               % left hand side
y_rhs   = rhs( beta_range, a, b, n0, n1, k0 );        % right hand side

% DEBUG plot dispersion
figure;
plot( beta_range, real(y_lhs) ); hold on;
plot( beta_range, real(y_rhs) );
legend('left hand side', 'right hand side');
xlabel('\beta (rad/nm)'); ylabel('dispersion');
title('DEBUG dispersion relation');
makeFigureNice();

% DEBUG plot dispersion, imag
figure;
plot( beta_range, imag(y_lhs) ); hold on;
plot( beta_range, imag(y_rhs) );
legend('left hand side', 'right hand side');
xlabel('\beta (rad/nm)'); ylabel('dispersion');
title('DEBUG dispersion relation, imaginary');
makeFigureNice();

% DEBUG plot abs(error)
figure;
plot( beta_range, abs(y_lhs - y_rhs) );
xlabel('\beta (rad/nm)'); ylabel('dispersion');
title('DEBUG dispersion relation, abs(error)');
makeFigureNice();

% solve for intersection of lhs and rhs
xlims   = (2*pi/lambda0) .* [ n1, n0 ];
% beta    = fzero( @(x) abs( lhs( x, a, n0, k0 ) - rhs( x, a, b, n0, n1, k0 ) ), xlims );
beta    = fminbnd( @(x) abs( lhs( x, a, n0, k0 ) - rhs( x, a, b, n0, n1, k0 ) ), xlims(1), xlims(2) );

% DEBUG calc err
err = abs( lhs( beta, a, n0, k0 ) - rhs( beta, a, b, n0, n1, k0 ) )

end

% -------------------------------------------------------------------------
% aux functions
% -------------------------------------------------------------------------

% left hand side of dispersion equation
function y = lhs( beta, a, n0, k0 )
    % inputs:
    %   beta    - vector of beta to solve for (rad/nm)
    %   a       - width of core, nm
    %   n0      - index of core and outer cladding
    %   k0      - 2*pi/lambda0, rad/nm
    
    kx  = sqrt( (k0^2)*(n0^2) - beta.^2 );
    y   = tan( kx*a );   
end

% right hand side of dispersion equation
function y = rhs( beta, a, b, n0, n1, k0 )
    % inputs:
    %   beta    - vector of beta to solve for (rad/nm)
    %   a       - width of core, nm
    %   b       - width of inter cladding, nm
    %   n0      - index of core and outer cladding
    %   n1      - index of inter cladding
    %   k0      - 2*pi/lambda0, rad/nms
    
    kx      = sqrt( (k0^2)*(n0^2) - beta.^2 );
    alpha   = sqrt( beta.^2 - (k0^2)*(n1^2) );
    y       = (alpha./kx) .* tanh( atanh( -1i * kx./alpha ) + alpha.*(b - a) );  
end




















