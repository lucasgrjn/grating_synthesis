function [] = save_fig_multiformat( h_fig, path, name, save_on )
% saves a matlab figure to the following formats:
% .eps, .png, .fig
% inputs:
%   h_fig = handle to figure
%   path = path to save (does not need last filesep)
%   name = name of file (no need to enter extension)
%     save_on = boolean
%                 if false, turn saving off (saves effort of writing if statements)

if nargin < 4
    save_on = true;
end

if save_on == true
    name = [path, filesep, name];

    %save .fig
    saveas(h_fig, name, 'fig');

    %save .png
    saveas(h_fig, name, 'png');

    %save .eps color
    saveas(h_fig, name, 'epsc');
end

end