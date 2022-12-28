function [x,y] = f_make_ellipse( x_intersect0, start_grating_width, start_grating_pos, ...
                                 neff, nclad, theta, npoints )
% returns x y coordinates of ellipse
%
% inputs:
%   x_intersect0 is the x axis intersect
%   npoints
%       type: scalar, int
%       desc: number of points to use
%   theta in radians

% get ellipse parameters
cur_taper_width = x_intersect0 * start_grating_width/start_grating_pos;
x0  = x_intersect0 / ( 1 + neff./( nclad * sin(theta) ) );     % center of ellipse
a   = x0 * neff / (  nclad * sin(theta) );                      % major axis of ellipse
b   = x0 * sqrt( neff^2 - ( nclad * sin(theta) )^2 ) / (  nclad * sin(theta) ); % minor axis

% make ellipse
y = linspace( -cur_taper_width/2, cur_taper_width/2, npoints );  % is 100 datapoints enough?
x = a * sqrt( 1 - (y./b).^2 ) + x0;

end