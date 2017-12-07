function [N,x,y]=couplingmesh(Lx,Ly,l,g,w,m,dx,dy,nBody,nBkgrnd)
%Lx x xip dimension 
%Ly y xip dimension 
%l coupling length
%g gap coupling distance
%w waveguide width
%dx,dy discretization
%nBody refractive index waveguide
%nBkgrnd refractive index enviroment
%m upper margin distance
% example Lx=100; Ly=50; l=30; g=0.5; w = 0.3; m=4; dx= 0.01; dy = 0.01; 
ds=Ly-2*(w+m)+g;
ls=0.5*(Lx-l);
xx1 = 0; xx2 = ls; yy1 = 0; yy2 = ds/2-g;
xVec = (xx1+dx/2):dx:(xx2-dx/2 + dx/2);  % added dx/2 to end to make *sure* last pixel is placed, true value is: (x1+dx/2):dx:(x2-dx/2);
yVec = (yy1+dy/2):dy:(yy2+w-dy/2 + dy/2);  % same story as line above
[xx,yy]=meshgrid(xVec,yVec);
f = @(xVec) ((yy2-yy1)/2)*cos(pi*(xVec-xx1)/(xx2-xx1)+pi)+(yy1+yy2)/2;
ys=f(xVec) ;
ysp=ys+w./cos(atan(dif1(f,xVec)));
n=ones(size(xx))*nBkgrnd;
n(yy > ones(length(yVec),1)*ys & yy< ones(length(yVec),1)*(ysp))=nBody;

c=ones(g/dy,Lx/dx)*nBkgrnd;
m=ones(m/dy,Lx/dx)*nBkgrnd;
wg=ones(length(yVec),l/dx)*nBkgrnd;
wg(1:w/dy,1:l/dx)=nBody;

N=[m;flipud([flipud(n),wg,n]);c;flipud(n),wg,n;m];
x= dx/2:dx:Lx;
y= dy/2:dy:Ly;
end

function df=dif1(F,x)
        h=1e-7;
        f2=feval(F,x+h);f1=feval(F,x-h);
        df=(f2-f1)/(2*h);
end