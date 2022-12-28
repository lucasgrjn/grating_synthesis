function [ field_overlap ] = f_field_position_overlap_fullvec_2D( E1, H1, E2, H2, normal_dir )
% calculates mode overlap vs. position using xcorr
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
    overlaps_E1_H2 = ifftshift( ifft2( conj( fft2( fftshift( E1.y ) ) ) .* fft2( fftshift( H2.z ) ) ) ) - ...
                     ifftshift( ifft2( conj( fft2( fftshift( E1.z ) ) ) .* fft2( fftshift( H2.y ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft2( fft2( fftshift( E2.y ) ) .* conj( fft2( fftshift( H1.z ) ) ) ) ) - ...
                     ifftshift( ifft2( fft2( fftshift( E2.z ) ) .* conj( fft2( fftshift( H1.y ) ) ) ) );
                 
elseif strcmp( normal_dir, 'y' )

    % first normalize both fields
    norm_factor_1   = E1.z(:)' * H1.x(:) - E1.x(:)' * H1.z(:) + H1.x(:)' * E1.z(:) - H1.z(:)' * E1.x(:);
    norm_factor_2   = E2.z(:)' * H2.x(:) - E2.x(:)' * H2.z(:) + H2.x(:)' * E2.z(:) - H2.z(:)' * E2.x(:);

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft2( conj( fft2( fftshift( E1.z ) ) ) .* fft2( fftshift( H2.x ) ) ) ) - ...
                     ifftshift( ifft2( conj( fft2( fftshift( E1.x ) ) ) .* fft2( fftshift( H2.z ) ) ) );
     
    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft2( fft2( fftshift( E2.z ) ) .* conj( fft2( fftshift( H1.x ) ) ) ) ) - ...
                     ifftshift( ifft2( fft2( fftshift( E2.x ) ) .* conj( fft2( fftshift( H1.z ) ) ) ) );

elseif strcmp( normal_dir, 'z' )

    % first normalize both fields
    norm_factor_1   = E1.x(:)' * H1.y(:) - E1.y(:)' * H1.x(:) + H1.y(:)' * E1.x(:) - H1.x(:)' * E1.y(:);
    norm_factor_2   = E2.x(:)' * H2.y(:) - E2.y(:)' * H2.x(:) + H2.y(:)' * E2.x(:) - H2.x(:)' * E2.y(:);

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft2( conj( fft2( fftshift( E1.x ) ) ) .* fft2( fftshift( H2.y ) ) ) ) - ...
                     ifftshift( ifft2( conj( fft2( fftshift( E1.y ) ) ) .* fft2( fftshift( H2.x ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft2( fft2( fftshift( E2.x ) ) .* conj( fft2( fftshift( H1.y ) ) ) ) ) - ...
                     ifftshift( ifft2( fft2( fftshift( E2.y ) ) .* conj( fft2( fftshift( H1.x ) ) ) ) );
                       
end

% calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
field_overlap = abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ ( norm_factor_1 .* norm_factor_2 );

end

