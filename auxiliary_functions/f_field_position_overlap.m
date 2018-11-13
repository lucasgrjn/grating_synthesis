function [ field_overlap ] = f_field_position_overlap( E1_z, H1_x, E2_z, H2_x )
% calculates mode overlap vs. position using xcorr
% currently assumes normal direction is y, E field polarized along z
% field is located along horizontal line parallel to x axis
%
% authors: bohan zhang
%
% Inputs:
%   E1_z
%       type: double, vector (column vector)
%       desc: electric field 1, z component
%   H1_x
%       type: double, vector (column vector)
%       desc: magnetic field 1, x component
%   E2_z
%       type: double, vector (column vector)
%       desc: electric field 2, z component
%   H2_x
%       type: double, vector (column vector)
%       desc: magnetic field 2, x component
%
% Outputs:
%   field_overlap
%       type: double, vector
%       desc: field overlap vs. position of fields 1 and 2


% first normalize both fields
norm_factor_1   = H1_x' * E1_z;
norm_factor_2   = H2_x' * E2_z;

% % DEBUG plot both fields
% figure;
% plot( 1:length( E1_z ), E1_z/sqrt(abs(norm_factor_1)) );
% title('Field E1, normalized');
% makeFigureNice();
% figure;
% plot( 1:length( E2_z ), E2_z/sqrt(abs(norm_factor_2)) );
% title('Field E2, normalized');
% makeFigureNice();
% figure;
% plot( 1:length( H1_x ), H1_x/sqrt(abs(norm_factor_1)) );
% title('Field H1, normalized');
% makeFigureNice();
% figure;
% plot( 1:length( H2_x ), H2_x/sqrt(abs(norm_factor_2)) );
% title('Field H2, normalized');
% makeFigureNice();

% calculate E1 and H2 xcorr
overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );

% calculate H1 and E2 xcorr
overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* fft( fftshift( H1_x ) ) ) );

% % DEBUG plot E1 and E2 in k space
% fft_E1 = fft( fftshift( E1_z ) );
% fft_E2 = fft( fftshift( E2_z ) );
% fft_H1 = fft( fftshift( H1_x ) );
% fft_H2 = fft( fftshift( H2_x ) );
% figure;
% plot( 1:length(fft_E1), ifftshift( fft_E1 ) );
% title('E1 in k space'); makeFigureNice();
% figure;
% plot( 1:length(fft_E2), ifftshift( fft_E2 ) );
% title('E2 in k space'); makeFigureNice();
% figure;
% plot( 1:length(fft_H1), ifftshift( fft_H1 ) );
% title('H1 in k space'); makeFigureNice();
% figure;
% plot( 1:length(fft_H2), ifftshift( fft_H2 ) );
% title('H2 in k space'); makeFigureNice();

% % DEBUG overplot E1 and H2 in freq domain
% figure;
% plot( 1:length(fft_E1), ifftshift( fft_E1 )/sqrt(abs(norm_factor_1)), 1:length(fft_H2), ifftshift( fft_H2 )/sqrt(abs(norm_factor_2)) );
% legend( 'e1', 'h2' );
% title('Overplot of E1 and H2 in k space (real)'); makeFigureNice();
% figure;
% plot( 1:length(fft_E1), ifftshift( abs(fft_E1) )/sqrt(abs(norm_factor_1)), 1:length(fft_H2), ifftshift( abs(fft_H2)/sqrt(abs(norm_factor_2)) ) );
% legend( 'e1', 'h2' );
% title('Overplot of E1 and H2 in k space (abs)'); makeFigureNice();


% % DEBUG plot overlaps E1 H2 and H1 E2
% figure;
% plot( 1:length( overlaps_E1_H2 ), overlaps_E1_H2 );
% title('Overlaps E1 H2');
% makeFigureNice();
% figure;
% plot( 1:length( overlaps_H1_E2 ), overlaps_H1_E2 );
% title('Overlaps H1 E2');
% makeFigureNice();

% calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );


% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v1');
% makeFigureNice()

% % calculate overlap function, no reals in denom
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( (norm_factor_1) .* (norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v2, no reals in denom');
% makeFigureNice()

% % calculate overlap function, abs * abs in denom
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v4, abs in denom');
% makeFigureNice()

% % calculate overlap function, switching order of conjugations
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( fft( fftshift( E2_z ) ) .* conj( fft( fftshift( H1_x ) ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v3, switching order of conjugations');
% makeFigureNice()

% % calculate overlap function, with just E1 and H2 overlap
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v5, just E1 and H2 overlap');
% makeFigureNice()

% % calculate overlap function, with just E2 and H1 overlap
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* fft( fftshift( H1_x ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_H1_E2 ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v6, just E2 and H1 overlap');
% makeFigureNice()

% % calculate overlap function, switching BOTH orders of conjugations
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( conj( fft( fftshift( E1_z ) ) ) .* ( fft( fftshift( H2_x ) ) ) ) );
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* ( fft( fftshift( H1_x ) ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v7, switching BOTH orders of conjugations');
% makeFigureNice()

% % Calculate overlap function, conjugating before the xcorr
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( conj(H2_x) ) ) ) ) );
% 
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( conj(E2_z) ) ) ) .* fft( fftshift( H1_x ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 + overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v8, conjugating before xcorr');
% makeFigureNice()

% % Changing the sign of the addition to subtraction
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );
% 
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* fft( fftshift( H1_x ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( overlaps_E1_H2 - overlaps_H1_E2 ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v9, subtraction instead of addition');
% makeFigureNice()

% % Absolute values in addition
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );
% 
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* fft( fftshift( H1_x ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v10, abs value of addition');
% makeFigureNice()

% % NO fftshifting
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifft( fft( E1_z ) .* conj( fft( H2_x ) ) );
% 
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifft( conj( fft( E2_z ) ) .* fft( H1_x ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ abs( real(norm_factor_1) .* real(norm_factor_2) );
% 
% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v11, no fftshifts');
% makeFigureNice()

% % Absolute values in addition and denom
% % calculate E1 and H2 xcorr
% overlaps_E1_H2 = ifftshift( ifft( fft( fftshift( E1_z ) ) .* conj( fft( fftshift( H2_x ) ) ) ) );
% 
% % calculate H1 and E2 xcorr
% overlaps_H1_E2 = ifftshift( ifft( conj( fft( fftshift( E2_z ) ) ) .* fft( fftshift( H1_x ) ) ) );
% 
% % calculate overlap function
% field_overlap = (1/4) * abs( abs(overlaps_E1_H2) + abs(overlaps_H1_E2) ).^2 ./ ( abs(norm_factor_1) .* abs(norm_factor_2) );

% % DEBUG plot field overlap
% figure;
% plot( 1:length(field_overlap), field_overlap );
% xlabel('x'); ylabel('overlap');
% title('Field overlap vs position v12, abs value of addition and denom');
% makeFigureNice()

end

