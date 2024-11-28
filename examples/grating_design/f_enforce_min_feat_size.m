function is_ok = f_enforce_min_feat_size( period, fill )
% Returns true if datapoint satisfies min. feature size demands, otherwise
% returns false if it violates
%
% Inputs:
%   period
%       type: array
%       desc: period of unit cell
%   fill
%       type: array
%       desc: ratio of layer length/period

% set these according to your process
min_line = 200;
min_space = 200;
 
is_ok = (period .* fill) >= min_line;
is_ok = is_ok & ( (period .* (1-fill)) >= min_space);
            
end