% unidirectional grating

clear; clc; close all;
tic;
% all units nm

% ------------------------------------------------------------------------|
% Paths and Save data struct/parameters
% ------------------------------------------------------------------------|

% folder containing refractive index data
addpath('refractiveindexdata');

% choose whether to save data
save_on = true;

% save directory (make sure to pick the right one)
save_dir = 'C:\Users\bozh7002\Google Drive\research\popovic group\data\gratings';   % drive folder on plab desktop

% filename is year_month_day hour_minute_second
filename = [datestr(clock,'yyyy_mm_dd HH_MM_SS'), ' data.mat'];
filename = [save_dir, '\', filename];

all_data = struct();
all_data.script_name = mfilename('fullpath');       % get name/path of script that made this data
all_data.run_time = datestr(clock);
all_data.notes = 'Bloch period 680nm';  % description of this simulation
all_data.units = 'nm';                              % units

% print some save status stuff
if save_on
    fprintf('Saving is on\n');
else
    fprintf('Saving is off\n');
end

% ------------------------------------------------------------------------|
% Parameters
% ------------------------------------------------------------------------|

lambda0 = 1310;     % in nm
k0 = 2*pi/lambda0;

% indices
temp_kelvin     = 273;  % room temperature
[n_cSi, ~, ~]   = index_Si_fits(lambda0/1000, temp_kelvin);
n_pSi           = index_IBM12SOI45nm_fits(lambda0/1000, 'polySi');
n_SiN           = 2.0;
[n_SiO2, ~, ~]  = index_SiO2_fits(lambda0/1000, temp_kelvin);

% dimensions
Lambda  = 680;       % grating period
w_t     = 370;       % width top of grating
w_b     = 330;       % width bottom of grating
t       = 75;        % thickness of grating teeth and SiN layer
w_o     = 150;       % offset width
t_air   = 700;       % thickness of bottom air layer
t_SiO2_bot  = 145;   % thickness of SiO2 layers bottom
t_SiO2_top  = 500;   % thickness of SiO2 layer top

% discretization 
dxz = 5;

% ------------------------------------------------------------------------|
% Make grating index (test structure)
% ------------------------------------------------------------------------|

[ n, x_vec, z_vec ] = make_grating_cell( Lambda, dxz, [ n_cSi, n_pSi, n_SiN, n_SiO2 ],...
                                    w_t, w_b, w_o, t, t_air, t_SiO2_bot, t_SiO2_top );

% display
figure; imagesc(z_vec, x_vec, n); set(gca,'Ydir','normal'); colorbar; axis equal;
colormap(jet); title('Index of unit cell');
xlabel('z (nm)'); ylabel('x (nm)');

% ------------------------------------------------------------------------|
% Run modesolver
% ------------------------------------------------------------------------|

% inputs
modes = 1;
BC = 0;         % 0 is pec, 1 is pmc

% pml options
h_pml           = 100;      % nm
stren_pml       = 600;     % strength of pml
poly_order_pml  = 5;        % polynomial order of pml
PML_options = [1, h_pml, stren_pml, poly_order_pml];

% sweeping guessk's
% n_guessk = 300;
guessk = linspace(0, 2*pi/Lambda, 50);
guessk = 7.1166e-04 + 3.6192e-04i;
% guessk = k0*0.9;

% sweeping offset
% w_o_all = linspace(0, 2*layer_offset, 50);
w_o = 150;

% runs will hold the data for each simulation
% annoying but i have to declare every single field within this struct
% before using it
all_data.runs = struct( 'lambda0', {}, 'width_top', {}, 'width_bot', [],...
                        'width_offset', {}, 'grating_thick', {}, 'thick_SiO2_bot', {},...
                        'thick_SiO2_top', {}, 'thick_air', {}, 'bloch_period', {},...
                        'dx', {}, 'dz', {}, 'x_vec', {}, 'z_vec', {},...
                        'n_modes', {}, 'BC', {}, 'k0', {}, 'guessk', {},...
                        'pml_opt', {}, 'field', {}, 'k', {}, 'n', {} );



i_run = 0;
n_runs = length(guessk)*length(w_o);

% Uncomment and run this loop for doing a sweep of pre-chosen guessk's and
% offsets
for i_guessk = 1:length(guessk)
    for i_w_o = 1:length(w_o)

        cur_run = struct();

        i_run = i_run + 1;
        fprintf('loop %i of %i\n',i_run,n_runs);

        % Make grating index
        [ n, x_vec, z_vec ] = make_grating_cell( Lambda, dxz, [ n_cSi, n_pSi, n_SiN, n_SiO2 ],...
                                    w_t, w_b, w_o(i_w_o), t, t_air, t_SiO2_bot, t_SiO2_top );

        % run modesolver
        try
            [Phi_1D, k] = complexk_mode_solver_2D_PML(n,dxz,k0,modes,guessk(i_guessk),BC,PML_options);

            % reshape field so dimensions are (x, z, mode#)
            field = reshape(Phi_1D, [fliplr(size(n)), size(Phi_1D,2)] );
            field = permute(field, [2 1 3]); % taking transpose of first two dimensions only
        catch
            % no modes found for this guessk
            Phi_1D = [];
            field = [];
            k = [];
        end


        cur_run.lambda0         = lambda0;
        cur_run.width_top       = w_t;
        cur_run.width_bot       = w_b;
        cur_run.width_offset    = w_o(i_w_o);
        cur_run.grating_thick   = t;
        cur_run.thick_SiO2_bot  = t_SiO2_bot;
        cur_run.thick_SiO2_top  = t_SiO2_top;
        cur_run.thick_air       = t_air;
        cur_run.bloch_period    = Lambda;
        cur_run.dx = dxz; cur_run.dz = dxz;
        cur_run.x_vec           = x_vec;
        cur_run.z_vec           = z_vec;
        cur_run.n_modes         = modes;
        cur_run.BC              = BC; % rewrite as string for legibility?
        cur_run.k0              = k0;
        cur_run.guessk          = guessk(i_guessk);
        cur_run.pml_opt         = PML_options;
        cur_run.field           = field;
        cur_run.k               = k;
        cur_run.n               = n;

        all_data.runs(end+1) = cur_run;
        
    end     % end offset loop
end     % end guessk loop




% save data file
if save_on
   fprintf('Saving...\n');
   save(filename,'all_data');
   fprintf('Saved.\n');
end

toc






