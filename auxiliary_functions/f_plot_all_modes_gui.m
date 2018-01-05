function [] = f_plot_all_modes_gui( phi, x, y, k )
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
%   x
%       type: double, vector
%       desc: x coordinates
%   y 
%       type: double, vector
%       desc: y coordinates
%   k
%       type: double, vector
%       desc: prop. eigenvalues vs. mode #

% Create a figure and axes
f   = figure('Visible','off');
ax  = axes('Units','pixels');

% default settings
mode_num    = 1;
mode_comp   = 'real';       % 'real', 'imag', 'amplitude', 'intensity'
phase_on    = false;        % determines whether the propagation phase is included

% two copies of field, one with and one without phase
phi_nophase = phi;
phi_phase   = phi;
for ii = 1:length(k)
    % for each mode, add in the phase
    phase_only          = repmat( exp( 1i * k(ii) * x ), length(y), 1 );
    phi_phase(:,:,ii)   = phi_phase(:,:,ii) .* phase_only;
end


% plot the first mode by default
plot_field();

% move axes over
ax.Position(1) = ax.Position(1) + 100;

% increase figure window size
f.Position(4) = f.Position(4) + 25;
f.Position(3) = f.Position(3) + 100;

% create drop down menu labels
nmodes          = size(phi,3);          % number of simulated modes
labelstrings    = {};                   % cell array of labels
for ii = 1:nmodes
    labelstrings{end+1} = [ ' mode ', num2str(ii) ];
end

% Create pop-up menu for selecting the mode
popup_selectmode = uicontrol('Style', 'popup',...
                           'String', labelstrings,...
                           'Position', [20, ax.Position(4), 100, 50],...
                           'Callback', @selectmode);   
                       
% Create pop-up menu for selecting component to plot
popup_selectcomponent = uicontrol('Style', 'popup',...
                                   'String', {' real', ' imag', ' amplitude', ' intensity'},...
                                   'Position', [20, ax.Position(4) - 25, 100, 50],...
                                   'Callback', @select_component);   
                               
% Create checkbox for selecting/deselecting phase
checkbox_phase = uicontrol('Style', 'checkbox',...
                                   'String', {'Phase on'},...
                                   'Position', [20, ax.Position(4) - 50, 100, 50],...
                                   'Callback', @toggle_phase);   

% Make figure visble after adding all components
f.Visible = 'on';

% function for selecting which mode to draw
function selectmode( source, event )
    % selects which mode to draw
    
    mode_num     = source.Value;

    plot_field();
end

% function for selecting which component to draw
function select_component( source, event )
    % selects which mode to draw
    
    mode_comp    = source.String{ source.Value };
    mode_comp    = mode_comp(2:end);                                        % remove whitespace prefix 

    plot_field();
end

% function for toggling phase on/off
function toggle_phase( source, event )

    if source.Value == true
        phi = phi_phase;
    else
        phi = phi_nophase;
    end 
    plot_field()
end

function plot_field()
    % depends on the variables:
    %   mode_num
    %   mode_comp
    % which are automatically updated by the other ui functions

    switch mode_comp
       
        case 'real'
            imagesc( x, y, real(phi(:,:,mode_num)) );
            title( sprintf( 'Mode %i, real', mode_num ));
            
        case 'imag'
            imagesc( x, y, imag(phi(:,:,mode_num)) );
            title( sprintf( 'Mode %i, imag', mode_num ));
            
        case 'amplitude'
            imagesc( x, y, abs(phi(:,:,mode_num)) );
            title( sprintf( 'Mode %i, amplitude', mode_num ));
            
        case 'intensity'
            imagesc( x, y, abs(phi(:,:,mode_num)).^2 );
            title( sprintf( 'Mode %i, intensity', mode_num ));
        
    end
    set( gca, 'ydir', 'normal' );
    colormap('redbluehilight');
    xlabel('x'); ylabel('y');
    colorbar;
    
end


end