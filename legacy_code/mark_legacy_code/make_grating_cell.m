function [ n, x_vec, z_vec ] = make_grating_cell( bloch_period, dxz, n_materials,...
                                    w_t, w_b, w_o, t, t_air, t_SiO2_bot, t_SiO2_top )
% Written by: bohan zhang
% 
% Description:
% makes a grating unit cell based on geometry in Mark's OI 2015 paper
% ALL units are nanometers
% 
% inputs:
%     bloch_period
%         : period of cell in Z (consequently is also the size of the Z domain)
%     dxz
%         : discretization (both x and z)
%     n_materials
%         : 1x4 array of material indices
%           [ n_cSi, n_pSi, n_SiN, n_SiO2 ]
%     w_t
%         : width of top tooth (in z)
%     w_b
%         : width of bottom tooth (in z)
%     w_o
%         : offset
%     t
%         : thickness of teeth and SiN layer (in x)
%     t_air
%         : thickness of bottom air layer
%     t_SiO2_bot
%         : thickness of SiO2 layer, bottom
%     t_SiO2_top
%         : thickness of SiO2 layer, top
% 
% outputs:
%     n
%         : unit cell index matrix
%     x_vec
%         : x coordinate array
%     z_vec
%         : z coordinate array

% separate material indexes
n_cSi = n_materials(1);
n_pSi = n_materials(2);
n_SiN = n_materials(3);
n_SiO2 = n_materials(4);

% size of space
X = t_air + t_SiO2_bot + t_SiO2_top + 3*t;
Z = bloch_period;

% Get # data points
nX = (X/dxz); nZ = (Z/dxz); % number of data points
if abs(round(nX) - nX) > eps(nX)
   error('Nx must be an integer\n');    % error, discretization is bad
end
if abs(round(nZ) - nZ) > eps(nZ)
   error('Nz must be an integer\n');
end

% create x, z coordinate vectors
x_vec = dxz.*(0:nX-1);    % discretization of the unit cell ( for plotting )
z_vec = dxz.*(0:nZ-1);    % discretization of the unit cell ( for plotting )


% index refraction, background is air
n = ones(nX,nZ);     % dimensions (x, z)

% add SiO2 above air
n( x_vec >= t_air, : ) = n_SiO2;

% make bottom tooth
n( x_vec >= t_air + t_SiO2_bot & x_vec < t + t_air + t_SiO2_bot,...
   z_vec >= w_o & z_vec < w_o + w_b ) = n_cSi;

% make top tooth
n( x_vec >= t + t_air + t_SiO2_bot & x_vec < 2*t + t_air + t_SiO2_bot,...
   z_vec < w_t ) = n_pSi;

% make SiN to right of top tooth
n( x_vec >= t + t_air + t_SiO2_bot & x_vec < 2*t + t_air + t_SiO2_bot,...
   z_vec >= w_t ) = n_SiN;

% make SiN above top tooth
n( x_vec >= 2*t + t_air + t_SiO2_bot & x_vec < 3*t + t_air + t_SiO2_bot,...
   z_vec < w_t + t ) = n_SiN;

% make SiN on the right edge of the domain
n( x_vec >= 2*t + t_air + t_SiO2_bot & x_vec < 3*t + t_air + t_SiO2_bot,...
   z_vec >= Z - t ) = n_SiN;

end

