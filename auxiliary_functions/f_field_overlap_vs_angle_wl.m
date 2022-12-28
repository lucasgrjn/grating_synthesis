function [ coupling_results ] = f_field_overlap_vs_angle_wl( results, thetas )
% calculate overlaps
%
% Inputs:
%   results
%       type: struct
%       desc: results, must include;
%               x, y, Ex, Ey, Hx, Hy, lambda, nclad, w0 (=MFD/2) (units m)
%
% Outputs:
%   coupling_results
%       type: struct
%       desc: holds these fields
%           field_overlap - overlap vs x vs y vs freq vs theta

% universal gaussian settings, hardcoding some things
d0              = 0;
nclad           = results.nclad;
w0              = results.w0;

y_fiber         = results.y - results.y(round(end/2));
x_fiber         = results.x - results.x(round(end/2));


field_overlap               = zeros( [ size(results.Ex), length(thetas) ] );            % dimensions x vs. y vs. freq vs. theta
max_overlap                 = zeros( size(results.Ex, 3), length(thetas) );             % dimensions max overlap vs. freq vs. theta
max_field_overlap           = zeros(size(results.Ex,1), size(results.Ex,2));            % dimensions x vs y
best_best_overlap           = 0;
pos_y_max_overlap_vs_theta  = zeros(size(thetas));                                      % best overlap position per theta

tic;
for i_theta = 1:length(thetas)
    % for each angle
    
    fprintf('coupling angle loop %i of %i\n', i_theta, length(thetas) );
    
    for ii = 1:size( results.Ex, 3 )
        % for each frequency (theres only one...)

%         fprintf('coupling loop %i of %i\n', ii, size(Ex,3) );

        % make fiber mode
        [E_fiber, H_fiber]  = f_fiberModeGaussian_2D( w0, results.lambda(ii), x_fiber, y_fiber, thetas(i_theta), d0, nclad );
       
        % assumign radiation upwards and TE polarized light
        E_fiber.y           = E_fiber.z;
        H_fiber.y           = H_fiber.z;
        
        field_overlap(:,:,ii,i_theta)     = f_field_position_overlap_fullvec_2D( ...
                                                    struct( 'x', results.Ex(:,:,ii), 'y', results.Ey(:,:,ii) ), ...
                                                    struct( 'x', results.Hx(:,:,ii), 'y', results.Hy(:,:,ii) ), ...
                                                    E_fiber, H_fiber, ...
                                                    'z');
        max_overlap(ii,i_theta)         = max( max( abs(field_overlap(:,:,ii,i_theta)) ) );

        toc;

    end     % end frequency loop

%     % get longitudinal position of max coupling
%     [~,indx_pos]                            = max( abs( field_overlap( round(end/2), : ) ) );
%     pos_y_max_overlap_vs_theta(i_theta)     = results.y(indx_pos);
    
    % save best results
    if max(max_overlap(:,i_theta)) > best_best_overlap
       
        best_theta          = thetas(i_theta);
        max_field_overlap   = field_overlap;
        best_best_overlap   = max_overlap(:,i_theta);
        
    end

end     % end theta loop

coupling_results = struct( 'field_overlap', field_overlap, ...
                            'max_overlap', max_overlap, ...
                            'pos_y_max_overlap_vs_theta', pos_y_max_overlap_vs_theta, ...
                            'best_theta', best_theta, ...
                            'max_field_overlap', max_field_overlap, ...
                            'best_best_overlap', best_best_overlap );

end

