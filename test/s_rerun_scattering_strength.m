% authors: bohan zhang
%
% temporary script for re-running scattering strength with improved
% calculation

clear; close all;

% load data
data = load( 'C:\Users\bz\Google Drive\research\popovic group\projects\grating synthesis\data\2018 01 15 grating synthesis data\obj.mat' );

% unpack data
v2struct(data);
v2struct(obj_as_struct);

% for each grating coupler, re-sim
i_loop = 0;
for ii = 1:size(GC_vs_fills,1)
    for jj = 1:size(GC_vs_fills,2)
        
        i_loop = i_loop + 1;
        fprintf('loop %i of %i\n', i_loop, numel(GC_vs_fills) );
        
        % grab grating cell
        GC = GC_vs_fills{ii,jj};
        
        % run sim
        GC = GC.runSimulation( 1, 0, [1, 200, 20, 2], GC.k );
        
        % re-save 
        GC_vs_fills{ii,jj}          = GC;
        scatter_str_vs_fills(ii,jj) = GC.alpha_down;    % assuming couplign direction was down
        
    end
end

% copy pasted from synthGrating
fill_tops       = fliplr( 0.3:0.025:0.95 );
fill_bots       = fliplr( 0.3:0.025:0.95 );

% re plot scattering strength alpha vs. fill
figure;
imagesc( fill_bots, fill_tops, real(scatter_str_vs_fills) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Scattering strength (real) vs. fill factors');

obj_as_struct.GC_vs_fills_v2            = GC_vs_fills;
obj_as_struct.scatter_str_vs_fills_v2   = scatter_str_vs_fills;

% scattering strength in 1/um
scatter_str_vs_fills_um = scatter_str_vs_fills * 1000;
figure;
imagesc( fill_bots, fill_tops, (scatter_str_vs_fills_um) );
colorbar; set( gca, 'ydir', 'normal' );
xlabel('bottom fill factor'); ylabel('top fill factor');
title('Scattering strength (1/um) vs. fill factors');