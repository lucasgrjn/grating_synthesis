function [E,H] = f_shape_gaussian_direction_pol( E, H, coupling_direction, pol )
% reshapes gaussian mode based on direction and polarization

switch pol

    case 'TE'

        switch coupling_direction
            case 'up'
                % radiation upwards
                E.y = E.z;
                H.y = H.z;
            case 'down'
                % radiation downwards
                E.y = E.z;
                H.y = -H.z;
        end

    case 'TM'
        
        switch coupling_direction
            case 'up'
                % radiation upwards
                tempex  = E.x;
                temphx  = H.x;
                E.x     = E.z;
                H.x     = H.z;
                E.y     = -tempex;
                H.y     = -temphx;
            case 'down'
                % radiation downwards
                tempex  = E.x;
                temphx  = H.x;
                E.x     = E.z;
                H.x     = H.z;
                E.y     = tempex;
                H.y     = temphx;
        end

end

end

