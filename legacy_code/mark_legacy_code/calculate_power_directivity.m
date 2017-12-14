function [ d ] = calculate_power_directivity( field, dxz, k, n_top_row, n_bot_row )
% calculating power down-up directivity
% 
% TODO:
%     add Sz into the mix too? - done
%     fix the derivatives so they are on the right grid. 
%     I am thinking do one forward, one backward derivative and take the average to put them back on the central grid
% 
% inputs:
%     field
%         : field distribution of unit cell, without phase factor (x,z)
%     dxz
%         : discretization in x and z
%     k
%         : e^jkz
%     n_top_row
%         : number of row to calculate power up
%     n_bot_row
%         : number of row to calculate power down

% calculate derivative of field
dx_field_f      = (field(3:end, :) - field(2:end-1, :))./dxz;   % forward derivative
dx_field_b      = (field(2:end-1, :) - field(1:end-2, :))./dxz; % backwards derivative
dx_field        = (dx_field_f + dx_field_b)./2;    % take average derivative
dx_field_conj   = conj(dx_field);                  

dz_field_f      = (field(:, 3:end) - field(:,2:end-1))./dxz;     % forward deriv
dz_field_b      = (field(:, 2:end-1) - field(:,1:end-2))./dxz;   % back deriv
dz_field        = (dz_field_f + dz_field_b)./2;    % average
dz_field_conj   = conj(dz_field);

% should I take half step of field?
% field_halfstep_x = (field(2:end, :) + field(1:end-1, :))/2;
% field_halfstep_z = (field(:, 2:end) + field(:, 1:end-1))/2;

% % DEBUGGING
% % plot dx field
% figure; set(gcf,'Position',[50 50 1700 900]);
% subplot(1, 2, 1);
% imagesc(imag(dx_field) ); set(gca, 'Ydir','normal');
% title('dx field averaged');
% % title( sprintf('dx Bloch envelope (imag), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
% % xlabel('z (nm)'); ylabel('x (nm)'); 
% colorbar; colormap(jet);
% subplot(1, 2, 2);
% imagesc(imag(dx_field_f) ); set(gca, 'Ydir','normal');
% title('dx field forward');

% % 
% subplot(1, 3, 2);
% imagesc(real(dx_field) ); set(gca, 'Ydir','normal'); 
% % title( sprintf('dx Bloch envelope (real), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
% % xlabel('z (nm)'); ylabel('x (nm)'); 
% colorbar; colormap(jet);
% 
% subplot(1, 3, 3);
% imagesc(abs(dx_field) ); set(gca, 'Ydir','normal'); 
% % title( sprintf('dx Bloch envelope (amplitude), mode %i of %i, k/b = %f', i_mode, length(k), k(i_mode)/b) );
% % xlabel('z (nm)'); ylabel('x (nm)'); 
% colorbar; colormap(jet);

% Sx_up = real(1i*field( n_top_row, : ).*dx_field_conj( n_top_row, : ));
% Sx_down = real(1i*field( n_bot_row, : ).*dx_field_conj(  n_bot_row, : ));

% trying halfstep avg thing, doesn't really change anything but makes
% indexing nicer
% Sx_up_half = real(1i*field_halfstep_x( n_top_row, : ).*dx_field_conj( n_top_row, : ));
% Sx_down_half = real(1i*field_halfstep_x( n_bot_row, : ).*dx_field_conj(  n_bot_row, : ));

% x component of poynting vector (real)
Sx_up = real(1i*field( n_top_row, : ).*dx_field_conj( n_top_row - 1, : ));
Sx_down = real(1i*field( n_bot_row, : ).*dx_field_conj(  n_bot_row - 1, : ));

% z component of poynting vector (real)
Sz_up = real( k*field( n_top_row, 2:end-1 ).*conj( field( n_top_row, 2:end-1 ) ) + ...
            1i*field( n_top_row, 2:end-1).*dz_field_conj( n_top_row, : ) );
Sz_down = real( k*field( n_bot_row, 2:end-1 ).*conj( field( n_bot_row, 2:end-1 ) ) + ...
            1i*field( n_bot_row, 2:end-1 ).*dz_field_conj( n_bot_row, :));

% magnitude of poynting vectors
S_up_mag = sqrt( Sx_up( 2:end-1).^2 + Sz_up.^2 );
S_down_mag = sqrt( Sx_down( 2:end-1 ).^2 + Sz_down.^2 );

P_up_mag = sum(S_up_mag(:));
P_down_mag = sum(S_down_mag(:));
        
% P_up_x = sum(Sx_up(:));
% P_down_x = sum(Sx_down(:));
% 
% P_up_z = sum(Sz_up(:));
% P_down_z = sum(Sz_down(:));
% 
% P_up_abs = abs(P_up_x) + abs(P_up_z);
% P_down_abs = abs(P_down_x) + abs(P_down_z);

d = P_down_mag/P_up_mag;

end












