function is_ok = f_enforce_min_feat_size( period, top_fill, bot_fill )
% Returns true if datapoint satisfies min. feature size demands, otherwise
% returns false if it violates
%
% User can add arguments that determine the min. feature sizes
%
% Inputs:
%   period
%       type: double, scalar
%       desc: period of unit cell
%   offset
%       type: double, scalar
%       desc: offset of bottom layer from left edge of cell
%   top_fill
%       type: double, scalar
%       desc: ratio of top layer length/period
%   bot_fill
%       type: double, scalar
%       desc: ratio of bot layer length/period

is_ok = true;

% hardcoded min. feature sizes (nm)
min_gap = 80;
min_wg  = 80;

% check for min. waveguide length is satisfied
top_len = top_fill * period;
bot_len = bot_fill * period;
if top_len < min_wg || bot_len < min_wg
    is_ok = false;
end

% check for min. gap is satisfied
top_gap = period - top_len;
bot_gap = period - bot_len;
if top_gap < min_gap || bot_gap < min_gap
    is_ok = false;
end
            
end