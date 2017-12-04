function [ N ] = f_indexBuilder( nlyrs, dlyrsx, dlyrsy, dx, dy)
% returns an index matrix (no averaging)
% Inputs: 
%   nlyrs:  m x n matrix describing index of refraction
%   dlyrsx: dimensions along columns of nlyrs
%   dlyrsy: dimensions along rows of nlyrs
%   dx:     x discretizaiton
%   dy:     y discretization
% Outputs: 
%   N: index matrix

% turn dimensions into pixels
dlyrsx_d = round(dlyrsx./dx);
dlyrsy_d = round(dlyrsy./dy);
N = zeros(sum(dlyrsy_d),sum(dlyrsx_d));
cumsum_x = [0 cumsum(dlyrsx_d)]; 
cumsum_y = [0 cumsum(dlyrsy_d)]; 

for II=1:size(nlyrs,1)
    for JJ=1:size(nlyrs,2)
        N(cumsum_y(II)+1:cumsum_y(II+1),cumsum_x(JJ)+1:cumsum_x(JJ+1)) = nlyrs(II,JJ);
    end    
end

end

