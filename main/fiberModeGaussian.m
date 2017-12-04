% Generate Gaussian-beam mode fields at a plane through z = d0 at angle of
% theta (in degrees)
%
%to be used for grating simulation overlaps
%
% (using H.A. Haus, Waves & Fields, Chapter 5)
% CMG November 21st, 2014
%
%
% Inputs:  w0  - 1/e beam RADIUS at waist
%          k0  - 2*pi/lambda(um)
%          pol - 0=Ex, 1=Ey
%          theta - angle in degrees
%          d0  - distance from the waist to intersect 
%                   - minus is focusing 
%                   - defaults to 0 (at waist) if no input
%
% Outputs: F.{Ex,Ey,Ez,Hx,Hy,Hz}

function [F] = fiberModeGaussian(w0, k0, pol, yvec, zvec, theta, d0)

%Default d0 to 0
if nargin < 7
    d0 = 0;                %Distance from waist (negative is focusing)
end

%Constants
c = 299792458;     %speed of light [m/s]
mu0 = 4e-7*pi;     %permeability of free space [m kg/(s^2 A^2)]
k0 = k0*10^6;      %1/m
omega = k0*c;      %angular frequency [hz]
w0 = w0*10^-6;     %[meters] radius
d0 = d0*10^-6;     %[meters] offset
    
%Convert to radians
theta = pi/180*theta;

%Rotate to fiber reference frame
zprimevec = (zvec*sin(theta))*10^-6 + d0;                                       %zprime [meters]
yprimevec = zvec*cos(theta)*10^-6;                                              %yprime [meters]
xprimevec = yvec*10^-6;                                                         %xprime [meters]

%Create mesh of prime coordinates
[yprime, xprime] = meshgrid(yprimevec, xprimevec); 
[zprime, temp] = meshgrid(zprimevec, xprimevec);

b = k0*w0^2/2;                                                                  %b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)


u00 = 1j.*sqrt(k0*b/pi).*(1./(zprime+1j.*b)).*...
    exp(-1j.*k0.*(xprime.^2+yprime.^2)./(2*(zprime+1j.*b)));                    %Equation (5.2) in Haus [1/meters]


if pol == 0
    %Haus calls TE where most of E field is in the x direction which we call xprime. Therefore xprime is the same as y 
    %zprime is in the longitudinal direction of the fiber
    Exprime = -1i.*omega*u00.*exp(-1j.*k0*zprime);                              %[1/(s meters)]                           
    Ezprime = -(xprime./(zprime+1j.*b)).*Exprime;                               %[1/(s meters)]                            
    
    Hyprime = -1i.*k0/mu0*u00.*exp(-1j.*k0*zprime);                             %[s^2 A^2/(meters^2 kg)]
    Hzprime = -(yprime./(zprime+1j.*b)).*Hyprime;                               %[s^2 A^2/(meters^2 kg)]
 
    tempEx = -1*Ezprime*cos(theta);
    if length(yvec)>1
        F.Ex = (tempEx(1:end-1,:)+tempEx(2:end,:))/2.;                          %change where field is defined since it is half pixel shift in the y direction
    else
        F.Ex = tempEx;   
    end
    F.Ey = Exprime;                                                             %This is the correct pixel alignment
    
    tempEz = -1*Ezprime*sin(theta);
    if length(yvec)>1
        temp2Ez = (tempEz(:,1:end-1)+tempEz(:,2:end))/2.;
        F.Ez = (temp2Ez(1:end-1,:)+temp2Ez(2:end,:))/2.;                        %change where field is defined since it is half pixel shift in y and z direction
    else 
        F.Ez = tempEz;
    end
    
    tempHx = Hyprime*sin(theta)-Hzprime*cos(theta);
    if length(yvec)>1
        F.Hx = (tempHx(:,1:end-1)+tempHx(:,2:end))/2.;                          %change where field is defined since it is half pixel shift in z direction
    else 
        F.Hx = tempHx;
    end
    
    F.Hy = zeros(size(F.Ez)); 
    F.Hz = -Hyprime*cos(theta)-Hzprime*sin(theta);                              %Off by half pixel in the x direciton
    
elseif pol == 1
    
    %For TM: Exprime -> Eyprime
    %        Eyprime -> -1*Exprime  which is 0 
    %        Ezprime -> Ezprime 
    %        Hxprime -> Hyprime     which is 0 
    %        Hyprime -> -1*Hxprime
    %        Hzprime -> Hzprime 
    
    Eyprime = -1j.*omega*u00.*exp(-1j.*k0*zprime);
    Ezprime = -(xprime./(zprime+1j.*b)).*Eyprime;
    
    Hxprime = 1j.*k0/mu0*u00.*exp(-1j.*k0*zprime);
    Hzprime = -(yprime./(zprime+1j.*b)).*Hxprime;
    
    tempHx = -1*Hzprime*cos(theta);
    if length(yvec)>1
        F.Hx = (tempHx(1:end-1,:)+tempHx(2:end,:))/2.;
    else 
        F.Hx = tempHx;
    end
    
    F.Hy = Hxprime;                                                             %This is the correct pixel alignment
    
    tempHz = -1*Hzprime*sin(theta);
    if length(yvec)>1
    temp2Hz = (tempHz(:,1:end-1)+tempHz(:,2:end))/2.;
    F.Hz = (temp2Hz(1:end-1,:)+temp2Hz(2:end,:))/2.;
    else
        F.Hz = tempHz; 
    end
    
    tempEx = -1*Ezprime*cos(theta)+Eyprime*sin(theta);
    if length(yvec)>1
        F.Ex = (tempEx(:,1:end-1)+tempEx(:,2:end))/2.;
    else
        F.Ex = tempEx;
    end
    
    F.Ey = zeros(size(F.Hz));
    F.Ez = -1*Ezprime*sin(theta)-Eyprime*cos(theta);                            %Off by half pixel in the x direciton
    
else
    error('Invalid Polarization');
end

end

