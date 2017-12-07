% Generate Gaussian-beam mode fields with waist w0 at z = 0.
% (using H.A. Haus, Waves & Fields, Chap. 5)
% MP June 16, 2006; taken from gaussbeammode.m of Aug 30, 2005.
%
% [Ez,Hx,neff] = gaussbeam1d(w0, k, nn, x)
%
% Inputs:  w0  - 1/e beam radius at waist
%          k   - k = n*k0 medium k-vector
%          nn  - index of medium
% Outputs: Ez,Hx fields and effective index, neff.
%
% NOTE: E,H output fields normalized so that sum(1/2 real(Ez Hx*)) = 1

function [Ez,Hx,neff] = gaussbeam1d(w0, k, nn, xvec)
Z0 = 299792458*4e-7*pi;     % 377 ohms free space impedance

yvec=0;
u00 = sqrt(2/pi)/w0 * exp(-(xvec.^2 + yvec.^2)/w0^2);
neff = 1 - 2/(k*w0)^2;     % Incorrect for non-air?  Find effective index (below vacuum index), due to finite transverse mode extent.
%Ez = i*k * u00;
%Hx = neff/Z0 * Ez;

%A = (2*pi)^(1/4) * Z0/(k*sqrt(neff*w0));  % Makes 1/(2*Zo) * sum(Ez'*Hx)*dx = 1
%A = pi^(-1/4) * Z0/(k/neff)*sqrt(neff*w0));  % Makes 1/(2*Zo) * sum(Ez'*Hx)*dx = 1
A = (2*pi)^(1/4) /k *sqrt(Z0*w0/nn);  % Makes 1/(2*Zo) * sum(Ez'*Hx)*dx = 1

Ez = A* 1i*k * u00;
Hx = nn/Z0 * Ez;
