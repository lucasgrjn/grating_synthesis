function [ N ] = f_twoLevelBuilder( nbackground,domain, nlyrs, dims, midpoint, dx, dy)
% builds two level grating

N = nbackground*ones(round(domain(1)/dx), round(domain(2)/dy)).';

% discretize inputs
midpoints_d = [round(midpoint(:,1)./dx), round(midpoint(:,2)./dy)];
dims_d = [round(dims(:,1)./dx), round(dims(:,2)./dy)];

% place structures
for II=1:length(nlyrs)
    startIndX = round(midpoints_d(II,1)-dims_d(II,1)/2);  
    endIndX = round(midpoints_d(II,1)+dims_d(II,1)/2);
    startIndY = round(midpoints_d(II,2));
    endIndY = round(midpoints_d(II,2)+dims_d(II,2));
    N(startIndY:endIndY,startIndX:endIndX) = nlyrs(II);
end

end

