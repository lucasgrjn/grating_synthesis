function [ field_overlap ] = f_field_position_overlap_fullvec_2D_notxcorr( E1, H1, E2, H2, normal_dir )
% calculates mode overlap vs. position using brute force sliding
% this function doesn't seem to work
% Fully vectorial overlap of fields on 2D grid
%
% I'm pretty sure this works in 1D too
%
% authors: bohan zhang
%
% All vectors must be [ n x m ]
%
% also i used my own annoying convention where im looking sideways at the waveguide
% instead of into the direction of propagation
%
% Inputs:
%   E1
%       type: struct
%       desc: The first field in the overlap
%             Has fields: x, y, z, where each field is the corresponding component of
%             E vs. 1D grid coordinates
%   H1
%       type: struct
%       desc: The first field in the overlap
%             Has fields: x, y, z, where each field is the corresponding component of
%             H vs. 1D grid coordinates (same grid as E1)
%   E2
%       type: struct
%       desc: The second field in the overlap
%             Has fields: x, y, z, where each field is the corresponding component of
%             E vs. 1D grid coordinates (does not have to be on same grid as E1, i think)
%   H2
%       type: struct
%       desc: The second field in the overlap
%             Has fields: x, y, z, where each field is the corresponding component of
%             H vs. 1D grid coordinates (does not have to be on same grid as JH1, i think)
%   normal_dir
%       type: string
%       desc: direction of normal to cross section, either 'x', 'y', or 'z'
%
% Outputs:
%   field_overlap
%       type: double, matrix
%       desc: field overlap vs. position of fields 1 and 2


if strcmp( normal_dir, 'x' )

    % first normalize both fields
    norm_factor_1   = E1.y(:)' * H1.z(:) - E1.z(:)' * H1.y(:) + H1.z(:)' * E1.y(:) - H1.y(:)' * E1.z(:);
    norm_factor_2   = E2.y(:)' * H2.z(:) - E2.z(:)' * H2.y(:) + H2.z(:)' * E2.y(:) - H2.y(:)' * E2.z(:);

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft2( fft2( fftshift( conj( E1.y ) ) ) .* conj( fft2( fftshift( H2.z ) ) ) ) ) - ...
                     ifftshift( ifft2( fft2( fftshift( conj( E1.z ) ) ) .* conj( fft2( fftshift( H2.y ) ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft2( conj( fft2( fftshift( E2.y ) ) ) .* fft2( fftshift( conj(H1.z) ) ) ) ) - ...
                     ifftshift( ifft2( conj( fft2( fftshift( E2.z ) ) ) .* fft2( fftshift( conj(H1.y) ) ) ) );

%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft2( fft2( fftshift( ( E1.y ) ) ) .* conj( fft2( fftshift( H2.z ) ) ) ) ) - ...
%                      ifftshift( ifft2( fft2( fftshift( ( E1.z ) ) ) .* conj( fft2( fftshift( H2.y ) ) ) ) );
% 
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft2( conj( fft2( fftshift( E2.y ) ) ) .* fft2( fftshift( (H1.z) ) ) ) ) - ...
%                      ifftshift( ifft2( conj( fft2( fftshift( E2.z ) ) ) .* fft2( fftshift( (H1.y) ) ) ) );
    
elseif strcmp( normal_dir, 'y' )

    % first normalize both fields
    norm_factor_1   = E1.z(:)' * H1.x(:) - E1.x(:)' * H1.z(:) + H1.x(:)' * E1.z(:) - H1.z(:)' * E1.x(:);
    norm_factor_2   = E2.z(:)' * H2.x(:) - E2.x(:)' * H2.z(:) + H2.x(:)' * E2.z(:) - H2.z(:)' * E2.x(:);

%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft2( conj( fft2( fftshift( E1.z ) ) ) .* fft2( fftshift( H2.x ) ) ) ) - ...
%                      ifftshift( ifft2( conj( fft2( fftshift( E1.x ) ) ) .* fft2( fftshift( H2.z ) ) ) );
%                  
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft2( conj( fft2( fftshift( conj(E2.z) ) ) ) .* fft2( fftshift( conj(H1.x) ) ) ) ) - ...
%                      ifftshift( ifft2( conj( fft2( fftshift( conj(E2.x) ) ) ) .* fft2( fftshift( conj(H1.z) ) ) ) );
    
    % get original size of the field
    origsize        = size(E1.x);
    overlaps_E1_H2  = zeros(origsize);
    overlaps_H1_E2  = zeros(origsize);
                 
    % pad all components
    E1z_pad = padarray(E1.z, size(E1.z));
    E1x_pad = padarray(E1.x, size(E1.x));
    E2z_pad = padarray(E2.z, size(E2.z));
    E2x_pad = padarray(E2.x, size(E2.x));
    H1z_pad = padarray(H1.z, size(H1.z));
    H1x_pad = padarray(H1.x, size(H1.x));
    H2z_pad = padarray(H2.z, size(H2.z));
    H2x_pad = padarray(H2.x, size(H2.x));
    
    % loop and calculate overlaps
    for i_z = 1 : (size(E1z_pad,1) - origsize(1))
        for i_x = 1 : (size(E1z_pad,2) - origsize(2))
            
%             E1_z_shift_conj = circshift( conj( E1.z ), -origsize/2 + [ i_z-1, i_x-1 ] );
%             E1_x_shift_conj = circshift( conj( E1.x ), -origsize/2 + [ i_z-1, i_x-1 ] );
            overlaps_E1_H2(i_z, i_x) = sum( conj( E1.z ) .* H2x_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ), 'all' ) ...
                                        - sum( conj( E1.x ) .* H2z_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ), 'all' );
            
%             overlaps_H1_E2(i_z, i_x) = sum( E2.z .* conj( H1x_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ) ), 'all' ) ...
%                                         - sum( E2.x .* conj( H1z_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ) ), 'all' );
              
            overlaps_H1_E2(i_z, i_x) = sum( E2.z .* conj( H1x_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ) ), 'all' ) ...
                                        - sum( E2.x .* conj( H1z_pad( i_z:i_z + origsize(1) - 1, i_x:i_x + origsize(2) - 1 ) ), 'all' );
                                    
        end
    end
                 
    
elseif strcmp( normal_dir, 'z' )

    % first normalize both fields
    norm_factor_1   = E1.x(:)' * H1.y(:) - E1.y(:)' * H1.x(:) + H1.y(:)' * E1.x(:) - H1.x(:)' * E1.y(:);
    norm_factor_2   = E2.x(:)' * H2.y(:) - E2.y(:)' * H2.x(:) + H2.y(:)' * E2.x(:) - H2.x(:)' * E2.y(:);

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft2( fft2( fftshift( conj( E1.x ) ) ) .* conj( fft2( fftshift( H2.y ) ) ) ) ) - ...
                     ifftshift( ifft2( fft2( fftshift( conj( E1.y ) ) ) .* conj( fft2( fftshift( H2.x ) ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft2( conj( fft2( fftshift( E2.x ) ) ) .* fft2( fftshift( conj(H1.y) ) ) ) ) - ...
                     ifftshift( ifft2( conj( fft2( fftshift( E2.y ) ) ) .* fft2( fftshift( conj(H1.x) ) ) ) );
                 
%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft2( fft2( fftshift( ( E1.x ) ) ) .* conj( fft2( fftshift( H2.y ) ) ) ) ) - ...
%                      ifftshift( ifft2( fft2( fftshift( ( E1.y ) ) ) .* conj( fft2( fftshift( H2.x ) ) ) ) );
% 
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft2( conj( fft2( fftshift( E2.x ) ) ) .* fft2( fftshift( (H1.y) ) ) ) ) - ...
%                      ifftshift( ifft2( conj( fft2( fftshift( E2.y ) ) ) .* fft2( fftshift( (H1.x) ) ) ) );
                 
end

% calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
field_overlap = abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ ( norm_factor_1 .* norm_factor_2 );

end

