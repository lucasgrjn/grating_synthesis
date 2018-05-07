% drawing contours
%
% requires that you preload a synth object

close all;

% directivity vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, 10*log10(Q.directivities_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Directivity (dB) vs. fill factors');


% plot contour
figure;
contour( 10*log10(Q.directivities_vs_fills) );

% calculate gradient, no interp
[ gradx, grady ] = gradient( 10*log10(Q.directivities_vs_fills) );

% plot gradient as vector field
figure;
quiver( Q.fill_top_bot_ratio, Q.fill_bots, gradx, grady );

% plot magnitude of gradient
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, sqrt( gradx.^2 + grady.^2 ) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Gradient of Directivity (mag) vs. fill factors');

% plot contour of gradient
figure;
contour( sqrt( gradx.^2 + grady.^2 ) );
title('Contour of mag of gradient of directivity (dB)');


% try filtering/smoothing before taking gradient and contours
filt_window              = [ 0.2, 0.5, 0.2; ...
                             0.5, 2.0, 0.5; ...
                             0.2, 0.5, 0.2 ];
filt_window = filt_window/sum(filt_window(:));
dir_vs_fills_dB_filtered = imfilter( 10*log10(Q.directivities_vs_fills), filt_window );

% directivity, filtered, vs. fill
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, dir_vs_fills_dB_filtered );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Directivity (dB), filtered, vs. fill factors');

% calc gradient with interp
[ gradx, grady ] = gradient( dir_vs_fills_dB_filtered );

% plot magnitude of gradient
figure;
imagesc( Q.fill_top_bot_ratio, Q.fill_bots, sqrt( gradx.^2 + grady.^2 ) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('top/bottom fill ratio'); ylabel('bottom fill factor');
title('Gradient of Directivity, filtered, (mag) vs. fill factors');

% plot contour of gradient
figure;
contour( sqrt( gradx.^2 + grady.^2 ) );
title('Contour of mag of gradient of directivity, filtered (dB)');