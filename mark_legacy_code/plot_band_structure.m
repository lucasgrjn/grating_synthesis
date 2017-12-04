% Plots band structure
% todo: plot band structure
% plot band structure for each file

clear; clc; close all;

% ------------------------------------------------------------------------|
% setup
% ------------------------------------------------------------------------|

% pick datafiles to load
PathName = 'C:\Users\bozh7002\Google Drive\research\popovic group\data\gratings\';  %plab desktop
FileNames = {'2017_01_21 18_12_36 data pml 100 600 2', ...
            '2017_01_22 14_04_04 data pml 100 1000 2', ......
            '2017_01_22 18_16_30 data pml 100 600 5' };
        
% plot labeling
legendstrs = { '100nm, 600 str, 2 order', ...
                '100nm, 1000 str, 2 order', ...
                '100nm, 600 str, 5 order' };
plot_formats = { 'bo', 'rx', 'g+', 'm*' };

% ------------------------------------------------------------------------|
% Plotting
% ------------------------------------------------------------------------|

b = 2*pi/650; % hardcoded for now, future runs can load bloch period from the dataset
h_fig = figure; % handle to figure

for i_file = 1:length(FileNames)
    
    k_all = [];
    
    % load data
    load( [PathName,FileNames{i_file}] );

    n_runs = length(all_data.runs);
    
    for i_run = 1:n_runs

        % get current run data
        cur_data    = all_data.runs(i_run);
        k           = cur_data.k;

        k_all = [k_all; k];
        
    end % end runs
    
    % plotting
    figure(h_fig);
    plot( real(k_all)/b, imag(k_all), plot_formats{i_file} ); hold on;
    xlabel('real(k)/b'); ylabel('imag(k)'); title('Real vs. Imag k');
    grid on;
    
end

% add legend
figure(h_fig);
legend(legendstrs);