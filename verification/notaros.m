% Reproduction from Bohan Zhang original code
clear; close all;

% import code
addpath(['main']);         % main
addpath(['..' filesep 'main']);         % main

% initial settings
disc        = 10;
units       = 'nm';
lambda      = 1550;
index_clad  = 1.0;
domain      = [ 2000, 470 ];
numcells    = 10;


% make object
GC = c_twoLevelGratingCell(  'discretization',   disc, ...
                            'units',            units, ...
                            'lambda',           lambda, ...
                            'domain_size',      domain, ...
                            'background_index', index_clad, ...
                            'num_cells',        numcells );

% Add a layer
height_y    = 260;
min_y       = (domain(1)-height_y)/2;
index       = 3.4;
GC          = GC.addLayer( min_y, height_y, index );

% add first rectangle
width_x     = 100;
min_x       = 0;
min_y       = min_y+20;
height_y    = 240;
index       = index_clad;
GC          = GC.addRect( min_x, min_y, width_x, height_y, index );

% add second rectangle
width_x     = 150;
min_x       = 130;
min_y       = min_y + (240-60);
height_y    = 60;
index       = index_clad;
GC          = GC.addRect( min_x, min_y, width_x, height_y, index );

% DEBUG plot the index
GC.plotIndex();
