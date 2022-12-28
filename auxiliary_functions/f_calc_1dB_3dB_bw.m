function [lambda_1db_bw, lambda_3db_bw] = f_calc_1dB_3dB_bw(lambda, coupling_eff)
% calculates coupling efficiency bandwidths
%
% coupling eff. is in absolute units (0-1) not dB

% first interpolate the data
lambdainterp    = linspace(min(lambda), max(lambda), 1e4);
coupling_interp = interp1( lambda, coupling_eff, lambdainterp );

% position of peak coupling
[ peakcoup, indxpeak ] = max(coupling_interp);

% 1db bw

% lower bound
[ ~, indxlow ]  = min(abs( 10*log10(coupling_interp(1:indxpeak)) - (10*log10(peakcoup) - 1) ) );
lambda_low      = lambdainterp(indxlow);

% upper bound
[ ~, indxhi ]  = min(abs( 10*log10(coupling_interp(indxpeak:end)) - (10*log10(peakcoup) - 1) ) );
lambda_hi      = lambdainterp(indxpeak + indxhi - 1);

lambda_1db_bw  = lambda_hi - lambda_low;

% 3dB bw

% lower bound
[ ~, indxlow ]  = min(abs( 10*log10(coupling_interp(1:indxpeak)) - (10*log10(peakcoup) - 3) ) );
lambda_low      = lambdainterp(indxlow);

% upper bound
[ ~, indxhi ]  = min(abs( 10*log10(coupling_interp(indxpeak:end)) - (10*log10(peakcoup) - 3) ) );
lambda_hi      = lambdainterp(indxpeak + indxhi - 1);

lambda_3db_bw  = lambda_hi - lambda_low;

end

