% script for loading and combining parallel sweep data
% authors: bohan zhang
%
% currently loads all the .mat files in the data directory data_dir, then
% combines them by fill, 

clear; close all;

% data directory
data_dir = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave\2017 11 29 sweep fill';
% data_dir = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave\2017 11 28 sweep fill';
% data_dir = 'C:\Users\bz\Google Drive\research\popovic group\code\grating_synthesis\grating_synth_refactor\test_datasave\2017 11 20 sweep fill';
% data_dir = [ pwd filesep 'test_datasave' filesep 'sweep fill' ];

% get names of files in the directory
files       = dir( data_dir );
filenames   = extractfield(files, 'name');
data_files  = {};

% load all the mat files
data_all = {};
for ii = 1:length(filenames)
    % for each file
   
    % filename
    fname = filenames{ii};
    
    % load data if .mat file
    if contains( fname, '.mat' )
        data_all{end+1}     = load( [ data_dir, filesep, fname ] ); 
        data_files{end+1}   = fname;
    end
    
end


% DEBUG check which fills were ran
fills_all = [];
for ii = 1:length(data_all)
    
    sweep_obj = data_all{ii}.sweep_obj;
    cur_fills = sweep_obj.fill_vec;
    
    fills_all = [ fills_all, cur_fills ];
    
end

% sort fills
fills_all = sort( fills_all );

% plot fills
figure;
plot( 1:length(fills_all), fills_all, '-o' );
ylabel('fill'); title('DEBUG simulated fills');
makeFigureNice();


% combine all the data together
sweep_obj_new = data_all{1}.sweep_obj;  % init sweep obj
for ii = 2:length(data_all)
    
    cur_sweep_obj = data_all{ii}.sweep_obj;
    
    % combine the fill vecs
    sweep_obj_new.fill_vec = [ sweep_obj_new.fill_vec, cur_sweep_obj.fill_vec ];
    
    % combine the sweep results data
    % tensors have dimensions ( fill, ratio, period, offset )
    names = fieldnames( sweep_obj_new.sweep_results );
    for jj = 1:length(names)
        sweep_obj_new.sweep_results.(names{jj}) = cat(  1, ...
                                                    sweep_obj_new.sweep_results.(names{jj}), ...
                                                    cur_sweep_obj.sweep_results.(names{jj}) );
    end
    
end

% sort fills in ascending order
[ fill_sorted, indx_sort_fill ] = sort( sweep_obj_new.fill_vec );
sweep_obj_new.fill_vec          = fill_sorted;
% sort sweep results
names = fieldnames( sweep_obj_new.sweep_results ); 
for ii = 1:length(names)
    temp_result                             = sweep_obj_new.sweep_results.(names{ii});
    sweep_obj_new.sweep_results.(names{ii}) = temp_result( indx_sort_fill, :, :, : );
end

% % DEBUG show the sorted sweep results
% sweep_obj_new.sweep_results.fill_tensor(:, 1, 1, 1)
% sweep_obj_new.sweep_results.ratio_tensor(:, 1, 1, 1)

% update the documentation
sweep_obj_new.data_notes        = 'combined sweep fill';
sweep_obj_new.orig_data_files   = filenames;
sweep_obj_new.data_filename     = 'sweep fill combined';
sweep_obj                       = sweep_obj_new;

% save the combined data to new folder
new_data_folder = [ data_dir, filesep, 'sweep_fill' ];
mkdir( new_data_folder );
save([ new_data_folder, filesep, sweep_obj.data_filename ], 'sweep_obj');






























