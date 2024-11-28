% synthesizes bi-directional (single level) grating
%
% MAIN SCRIPT, the other bidir scripts are deprecated

clear; close all;

% ----------------------------
% dependencies
% grating synthesis codes
addpath(genpath('../../..'));

% ----------------------------
% main

% choose where to save data
loadsave_dir = pwd;

% inputs
lambda          = 1550;
optimal_angle   = 15;
disc            = 10;
apodized        = true; % select true to apodize, false for uniform
MFD = 10.4*1e3;
        
% name this grating design. each synthesized obj with design space gets a
% unique folder
gc_name = [ datestr(now,'yymmdd_HHMM') '_lambda' num2str(lambda) '_optangle' ...
            num2str(optimal_angle) '_dx' num2str(disc) '_MFD' num2str(MFD) ];
% make saving folder
loadsave_dir = [ loadsave_dir filesep gc_name ];
mkdir(loadsave_dir);

% first, generate the design space
% this will take some time
synth_obj = f_run_designspace_singlelevel( lambda, optimal_angle, disc );

% you can plot the design space to see what it looks like
figure('name', 'design_space');
yyaxis left;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, synth_obj.sweep_variables.scatter_str_vs_fill ); hold on;
ylabel('\alpha'); makeFigureNice();
yyaxis right;
plot( synth_obj.sweep_variables.fill_ratios_to_sweep, synth_obj.sweep_variables.periods_vs_fill );
ylabel('\Lambda'); makeFigureNice();
xlabel('DC');

% generate final design
input_wg_type = 'full'; % set this depending on what the input geometry looks like
                        % since we are creating the grating by etching, the
                        % input waveguide type is "full" e.g. starting from
                        % 100% duty cycle

% this function generates an apodized grating for fiber-coupling
synth_obj = synth_obj.generate_final_design_apodized_gaussian( MFD, input_wg_type, @f_enforce_min_feat_size );

% save the final synthesis object and workspace
save( [ loadsave_dir filesep 'synthobjfinal' ] );

% plot final design
f_plot_final_design( synth_obj, 'single' );






















