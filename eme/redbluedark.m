function J = redblue(m)
%REDBLUE    Color scheme with +ve values blue, -ve red, and zero white.
%   REDBLUE(M), a variant of JET(M), is an M-by-3 matrix containing
%   the default colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red.
%   JET, with no arguments, is the same length as the current colormap.
%   Use COLORMAP(JET).
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   August 8, 2003, Milos Popovic.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
n = floor(m/2);
if (mod(m,2) == 1) ODD = 1; else ODD = 0; end;
uL = [0:1:n-1]/(n-1 +ODD);
uR = fliplr(uL);
if (ODD) uC = 1; else uC = []; end;
uo = ones(1,n);

%r = [uo uC uR].';
%g = [uL uC uR].';
%b = [uL uC uo].';
r = [(uo+uL)/2 uC uR].';
g = [uL uC uR].';
b = [uL uC (uo+uR)/2].';
J = [r g b];

