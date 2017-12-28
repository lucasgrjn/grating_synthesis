function [] = f_plot_all_modes_gui( phi )
% authors: bohan zhang
%
% Function that plots all the modes from complex k modesolver, new ver, in
% a convenient gui
%
% Inputs:
%   phi
%       type: double, tensor
%       desc: Field (envelope), dimensions y vs x vs mode #, where y = in
%             plane transverse dimension and x = direction of propagation

% Create a figure and axes
f   = figure('Visible','off');
ax  = axes('Units','pixels');

% plot the first mode by default
imagesc( real(phi(:,:,1)) );
set( gca, 'ydir', 'normal' );
colorbar;
title('Mode 1, real');

% increase figure window size
f.Position(4) = f.Position(4) + 25;

% create drop down menu labels
nmodes          = size(phi,3);          % number of simulated modes
labelstrings    = {};                   % cell array of labels
for ii = 1:nmodes
    labelstrings{end+1} = [ 'mode ', num2str(ii) ];
end

% Create pop-up menu
popup = uicontrol('Style', 'popup',...
       'String', labelstrings,...
       'Position', [20, ax.Position(4) + 25, 100, 50],...
       'Callback', @selectmode);    

% Make figure visble after adding all components
f.Visible = 'on';

% function for selecting which mode to draw
function selectmode( source, event )
    % selects which mode to draw
    
    modenum     = source.Value;
    modename    = source.String;

    imagesc( real(phi(:,:,modenum)) );
    colorbar;
    title( sprintf( 'Mode %i, real', modenum ));
    
end


end