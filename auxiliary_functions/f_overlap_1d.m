function [max_overlap, pos_best] = f_overlap_1d(Ez, Hx, w0, lambda, theta, x, nclad )
% helper function, calculates overlap in 1D, picks out the wavelength that
% has the best overlap, and returns the overlap when the fiber is
% positioned for coupling at that wavelength
%
% only for coupling to gaussian
%
% look at D:\Google Drive\research\popovic group\papers\drafts\2020 07 20 -
% high eff grating\code\s_load_best_grating_2dfdtdresults.m for the inputs

d0              = 0;
unit_scaling    = 1.0;
x_fiber         = x - x(round(end/2));

% calculate overlap using function
field_overlap       = zeros( size(Ez) );                              % dimensions x vs. freq
max_overlap         = zeros( size(Ez, 2), 1 );
for ii = 1:size( Ez, 2 )
    
    % make fiber mode
    [Ez_fiber, Hx_fiber]    = f_fiberModeGaussian( w0, lambda(ii), x_fiber, unit_scaling, theta, d0, nclad );
    
    field_overlap(:,ii)     = f_field_position_overlap_fullvec_2D( ...
                                        struct( 'z', Ez(:,ii), 'x', zeros(size(Ez(:,ii))) ), ...
                                        struct( 'x', Hx(:,ii), 'z', zeros(size(Hx(:,ii))) ), ...
                                        struct( 'z', Ez_fiber, 'x', zeros(size(Ez_fiber)) ), ...
                                        struct( 'x', Hx_fiber, 'z', zeros(size(Hx_fiber)) ), ...
                                        'y');
    max_overlap(ii)         = max( abs(field_overlap(:,ii)) );
    
end

% pick out the coupling eff at the position for best coupling
[~, indx_wl_best]   = max(max_overlap);
[~, pos_best]       = max(field_overlap(:,indx_wl_best));
max_overlap         = field_overlap(pos_best,:);

end

