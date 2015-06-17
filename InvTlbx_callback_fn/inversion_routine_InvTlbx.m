% inversion routine that interacts with Christian's code

%% preparation

% number of iterations
InvProps.niter = 60;
istart = 1;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, ~, apply_hc, use_grav, fix_velocities, ...
    use_matfile_startingmodel, starting_model, bg_model_type,...
    true_model_type, ~, ~, ...
    parametrisation, param_plot, rec_g, ~, ~, ~, ...
    ~, ~, smoothgwid] = get_input_info;

% kernels to be added in same parametrisation as inversion.
param_addknls = parametrisation;

% project folder
output_path = ['./output/',project_name,'/'];
mkdir('./output/',project_name)

% save input_parameters to this path
savename = [output_path,project_name,'.input_parameters.m'];
copyfile('./input/input_parameters.m',savename)

%% welcome

disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp(['INVERSION RUN ', project_name]);
disp(['-- parametrisation:  ', parametrisation]);
disp '-- using seis? ..... YES!!'
if strcmp(use_grav , 'yes')
    disp '-- using grav? ..... YES!!'
else
    disp '-- using grav? ..... no'
end
if strcmp(apply_hc , 'yes')
    disp '-- using h.c.? ..... YES!!'
else
    disp '-- using h.c.? ..... no'
end
if strcmp(fix_velocities , 'yes')
    disp '-- fixing vels? .... YES!!'
else
    disp '-- fixing vels? .... no'
end
disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';


%% OBS

% freq consists of freqs.v_obs and freqs.frequency
if ~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var')
    if istart == 1
        disp 'no OBS present, preparing obs...';
        [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path,true_model_type);
    else
        error('iter > 1 but there are no observed properties!')
    end
else
    disp 'obs properties all present... proceeding...'
end


savename = [output_path,'/obs.all-vars.mat'];
if ~(exist(savename, 'file'))
    disp 'saving obs variables to matfile...'
    save(savename, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');
end

%=========================================================================

%% Model preparation

disp 'preparing inversion starting model...'

% save initial model (from input_parameters)
Model_start = update_model();

% save background model (for plotting purposes)
Model_bg = update_model(bg_model_type);

% set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));


% plot initial model
if istart == 1
    disp 'plotting initial model'
    fig_mod = plot_model_diff(Model_start,Model_bg,param_plot);
    set(fig_mod,'Renderer','painters');
    titel = [project_name,': starting model'];
    mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
    figname = [output_path,'/iter0.starting-model.',param_plot,'.png'];
    print(fig_mod,'-dpng','-r400',figname);
    close(fig_mod);
    
    
    % if 1st iter model @ matfile, load matfile
    if strcmp(use_matfile_startingmodel,'yes')
        load(starting_model)
        Model(1) = Model_out;
        clearvars Model_out;
    else
        Model(1) = Model_start;
    end
    
    
end

%% calculate initial model misfit per frequency

misfit_init = calc_initial_misfits(Model(1), sObsPerFreq, g_obs);
save([output_path,'/initial_misfits.mat'], 'misfit_init', '-v6');

%% set initial inversion values
whichFrq = 1;

sEventInfo = sObsPerFreq(whichFrq).sEventInfo;
sEventObs   = sObsPerFreq(whichFrq).sEventObs;


%% preparing input for the Inversion Toolbox

% perpare usr_par

usr_par.InvProps        = InvProps;
usr_par.Model(1)        = Model(1);
usr_par.sObsPerFreq     = sObsPerFreq;
usr_par.misfit_init     = misfit_init;
usr_par.whichFrq        = whichFrq;
usr_par.g_obs           = g_obs;
usr_par.sEventInfo      = sEventInfo;
usr_par.sEventObs       = sEventObs;
usr_par.Model_bg        = Model_bg;
usr_par.parametrisation = parametrisation;
usr_par.use_grav        = use_grav;
usr_par.smoothgwid      = smoothgwid;
usr_par.rec_g           = rec_g;
usr_par.output_path     = output_path;
if exist('Model_real', 'var')
    usr_par.Model_real  = Model_real;
end


%% calling onto the Inversion Toolbox

% call steepest gradient, L-BFGS, something else.

m = map_parameters_to_m(Model(1), usr_par);



