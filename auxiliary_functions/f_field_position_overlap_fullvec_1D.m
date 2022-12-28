function [ field_overlap, overlap_phase ] = f_field_position_overlap_fullvec_1D( E1, H1, E2, H2, normal_dir )
% calculates mode overlap vs. position using xcorr
% Fully vectorial overlap of fields on 1D grid
% I THINK OBSOLETE use 2d version instead
%
% authors: bohan zhang
%
% All vectors must be [ n x 1 ]
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
%       type: double, vector
%       desc: field overlap vs. position of fields 1 and 2

if strcmp( normal_dir, 'x' )
    
%     % pad arrays
%     [ E1.y, E2.y ] = pad_arrays( E1.y, E2.y );
%     [ E1.z, E2.z ] = pad_arrays( E1.z, E2.z );
%     [ H1.y, H2.y ] = pad_arrays( H1.y, H2.y );
%     [ H1.z, H2.z ] = pad_arrays( H1.z, H2.z );
    
    % first normalize both fields
    norm_factor_1   = H1.y' * E1.z;
    norm_factor_2   = H2.y' * E2.z;

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1.y ) ) .* conj( fft( fftshift( H2.z ) ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2.y ) ) ) .* fft( fftshift( H1.z ) ) ) );
    
elseif strcmp( normal_dir, 'y' )
    
%     % pad arrays
%     [ E1.x, E2.x ] = pad_arrays( E1.x, E2.x );
%     [ E1.z, E2.z ] = pad_arrays( E1.z, E2.z );
%     [ H1.x, H2.x ] = pad_arrays( H1.x, H2.x );
%     [ H1.z, H2.z ] = pad_arrays( H1.z, H2.z );

    % first normalize both fields
    norm_factor_1   = H1.x' * E1.z;
    norm_factor_2   = H2.x' * E2.z;
    
%     % first normalize both fields
%     norm_factor_1   = E1.z(:)' * H1.x(:) - E1.x(:)' * H1.z(:) + H1.x(:)' * E1.z(:) - H1.z(:)' * E1.x(:);
%     norm_factor_2   = E2.z(:)' * H2.x(:) - E2.x(:)' * H2.z(:) + H2.x(:)' * E2.z(:) - H2.z(:)' * E2.x(:);

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1.z ) ) .* conj( fft( fftshift( H2.x ) ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2.z ) ) ) .* fft( fftshift( H1.x ) ) ) );
    
elseif strcmp( normal_dir, 'z' )

%     % pad arrays
%     [ E1.y, E2.y ] = pad_arrays( E1.y, E2.y );
%     [ E1.z, E2.z ] = pad_arrays( E1.z, E2.z );
%     [ H1.y, H2.y ] = pad_arrays( H1.y, H2.y );
%     [ H1.z, H2.z ] = pad_arrays( H1.z, H2.z );
    
    % first normalize both fields
    norm_factor_1   = H1.x' * E1.y;
    norm_factor_2   = H2.x' * E2.y;

    % calculate E1 and H2 xcorr
    overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1.x ) ) .* conj( fft( fftshift( H2.y ) ) ) ) );

    % calculate H1 and E2 xcorr
    overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2.x ) ) ) .* fft( fftshift( H1.y ) ) ) );
    
end     % end if else statement

% calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
% field_overlap = abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ ( norm_factor_1 .* norm_factor_2 );
overlap_phase = angle( overlaps_E1_H2 + overlaps_H1_E2 ); % no idea if this is legit or not

% if strcmp( normal_dir, 'x' )
% 
%     % first normalize both fields
%     norm_factor_1   = E1.y(:)' * H1.z(:) - E1.z(:)' * H1.y(:) + H1.z(:)' * E1.y(:) - H1.y(:)' * E1.z(:);
%     norm_factor_2   = E2.y(:)' * H2.z(:) - E2.z(:)' * H2.y(:) + H2.z(:)' * E2.y(:) - H2.y(:)' * E2.z(:);
% 
%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( conj( E1.y ) ) ) .* conj( fft( fftshift( H2.z ) ) ) ) ) - ...
%                      ifftshift( ifft( fft( fftshift( conj( E1.z ) ) ) .* conj( fft( fftshift( H2.y ) ) ) ) );
% 
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2.y ) ) ) .* fft( fftshift( conj(H1.z) ) ) ) ) - ...
%                      ifftshift( ifft( conj( fft( fftshift( E2.z ) ) ) .* fft( fftshift( conj(H1.y) ) ) ) );
%     
% elseif strcmp( normal_dir, 'y' )
% 
%     % first normalize both fields
%     norm_factor_1   = E1.z(:)' * H1.x(:) - E1.x(:)' * H1.z(:) + H1.x(:)' * E1.z(:) - H1.z(:)' * E1.x(:);
%     norm_factor_2   = E2.z(:)' * H2.x(:) - E2.x(:)' * H2.z(:) + H2.x(:)' * E2.z(:) - H2.z(:)' * E2.x(:);
% 
%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( conj( E1.z ) ) ) .* conj( fft( fftshift( H2.x ) ) ) ) ) - ...
%                      ifftshift( ifft( fft( fftshift( conj( E1.x ) ) ) .* conj( fft( fftshift( H2.z ) ) ) ) );
%                  
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft( conj( fft2( fftshift( E2.z ) ) ) .* fft( fftshift( conj(H1.x) ) ) ) ) - ...
%                      ifftshift( ifft( conj( fft2( fftshift( E2.x ) ) ) .* fft( fftshift( conj(H1.z) ) ) ) );
% 
%     
% elseif strcmp( normal_dir, 'z' )
% 
%     % first normalize both fields
%     norm_factor_1   = E1.x(:)' * H1.y(:) - E1.y(:)' * H1.x(:) + H1.y(:)' * E1.x(:) - H1.x(:)' * E1.y(:);
%     norm_factor_2   = E2.x(:)' * H2.y(:) - E2.y(:)' * H2.x(:) + H2.y(:)' * E2.x(:) - H2.x(:)' * E2.y(:);
% 
%     % calculate E1 and H2 xcorr
%     overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( conj( E1.x ) ) ) .* conj( fft( fftshift( H2.y ) ) ) ) ) - ...
%                      ifftshift( ifft( fft( fftshift( conj( E1.y ) ) ) .* conj( fft( fftshift( H2.x ) ) ) ) );
% 
%     % calculate H1 and E2 xcorr
%     overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2.x ) ) ) .* fft( fftshift( conj(H1.y) ) ) ) ) - ...
%                      ifftshift( ifft( conj( fft( fftshift( E2.y ) ) ) .* fft( fftshift( conj(H1.x) ) ) ) );
%                  
% end
% 
% % calculate overlap function
% % field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% % field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
% field_overlap = abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ ( norm_factor_1 .* norm_factor_2 );

end     % end function


% aux functions

function [ array1, array2 ] = pad_arrays( array1, array2 )
% Zero pads two arrays of different length so they are the same length 
% so I can do FFTs properly
%
% expects inputs to be [ 1 x n ]

% get longer length of the two
longer_len = max( [ length(array1), length(array2) ] );

% pad array 1
pad_len     = longer_len - length(array1);
pad_low     = ceil( pad_len/2 );
pad_high    = floor( pad_len/2 );
array1      = [ zeros(1,pad_low), array1, zeros(1,pad_high) ];

% pad array 2
pad_len     = longer_len - length(array2);
pad_low     = ceil( pad_len/2 );
pad_high    = floor( pad_len/2 );
array2      = [ zeros(1,pad_low), array2, zeros(1,pad_high) ];

end

