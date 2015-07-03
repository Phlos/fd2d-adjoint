% inversion routine that interacts with Christian's code

%% preparation

% number of iterations
InvProps.niter = 60;
istart = 1;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, ~, apply_hc, use_grav, fix_velocities, ...
    use_matfile_startingmodel, starting_model, bg_model_type,...
    true_model_type, ~, change_freq_every, ...
    parametrisation, param_plot, rec_g, ~, ~, ~, ...
    ~, InvProps.stepInit, smoothgwid] = get_input_info;

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

obs_file = [output_path,'/obs.all-vars.mat'];
% freq consists of freqs.v_obs and freqs.frequency
if ((~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var'))  && exist(obs_file, 'file'))
    disp 'loading OBS file...'
    load(obs_file);
elseif (~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var'))
    if istart == 1
        disp 'no OBS present, preparing obs...';
        [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path,true_model_type);
        save(savename, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');
    else
        error('iter > 1 but there are no observed properties!')
    end
else
    disp 'obs properties all present... proceeding...'
end


% savename = [output_path,'/obs.all-vars.mat'];
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

% % set the background values for plot_model to mode of the initial model
% middle.rho    = mode(Model_start.rho(:));
% middle.mu     = mode(Model_start.mu(:));
% middle.lambda = mode(Model_start.lambda(:));


% plot initial model
if istart == 1
%     disp 'plotting initial model'
%     fig_mod = plot_model_diff(Model_start,Model_bg,param_plot);
%     set(fig_mod,'Renderer','painters');
%     titel = [project_name,': starting model'];
%     mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
%     figname = [output_path,'/iter0.starting-model.',param_plot,'.png'];
%     print(fig_mod,'-dpng','-r400',figname);
%     close(fig_mod);
    
    
    % if 1st iter model @ matfile, load matfile
    if strcmp(use_matfile_startingmodel,'yes')
        load(starting_model)
        Model(1) = Model_out;
        clearvars Model_out;
    else
        Model(1) = Model_start;
    end
    
%     % model
%     disp 'plotting initial model'
%     iter = 1;
%     fig_mod = plot_model_diff(Model(iter),Model_bg,param_plot);
%     titel = [project_name,': model diff of iter ', num2str(iter), ' and bg model'];
%     mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
%     figname = [output_path,'/iter',num2str(iter,'%03d'),'.model-diff.',param_plot,'.png'];
%     print(fig_mod,'-dpng','-r400',figname);
%     close(fig_mod);
%     clearvars('fig_mod');
    
    
end

%% calculate initial model misfit per frequency (if not present yet)

init_misfit_file = [output_path,'/initial_misfits.mat'];
% skip calculating initial misfits if they are already in the workspace
% or if there's a initial misfit file in the inversion directory
if  (~exist('misfit_init', 'var') && exist(init_misfit_file, 'file'))
    disp 'loading initial misfit file'
    load(init_misfit_file);
elseif ~exist('misfit_init', 'var')
    disp 'calculating initial misfits'
    misfit_init = calc_initial_misfits(Model(1), sObsPerFreq, g_obs);
    save(init_misfit_file, 'misfit_init', '-v6');
else
    disp 'initial misfits already present... proceeding...';
end

%% set initial inversion values

% set the stage for multiple frequency inv.
whichFrq = 1;
cfe = change_freq_every;
cumulative_iter = 1;
nfreq = numel(sObsPerFreq);

sEventInfo = sObsPerFreq(whichFrq).sEventInfo;
sEventObs   = sObsPerFreq(whichFrq).sEventObs;




%% preparing input for the Inversion Toolbox

% perpare usr_par

usr_par.InvProps        = InvProps;
usr_par.Model(1)        = Model(1);
usr_par.sObsPerFreq     = sObsPerFreq;
usr_par.misfit_init     = misfit_init;
% usr_par.whichFrq        = whichFrq;
usr_par.g_obs           = g_obs;
usr_par.sEventInfo      = sEventInfo;
usr_par.sEventObs       = sEventObs;
usr_par.Model_bg        = Model_bg;
usr_par.parametrisation = parametrisation;
usr_par.use_grav        = use_grav;
usr_par.smoothgwid      = smoothgwid;
usr_par.rec_g           = rec_g;
usr_par.output_path     = output_path;
usr_par.cumulative_iter = cumulative_iter;
if exist('Model_real', 'var')
    usr_par.Model_real  = Model_real;
end


%% calling onto the Inversion Toolbox

% remove old fwd fields
TempFolder = [output_path,'/fwd_temp'];
if (exist(TempFolder, 'dir'))
    rmdir([output_path,'/fwd_temp'], 's');
end

% call steepest gradient, L-BFGS, something else.

m = map_parameters_to_m(Model(1), usr_par);



output_log = [output_path,'/lbfgs_output_log.txt'];

% loop over L-BFGS calls at different frequencies
for ifreq = 1:nfreq
   
    disp '=============================';
    disp '  going to a new frequency   ';
    disp(['  (frequency ',num2str(ifreq),'/',num2str(nfreq)',')']);
    disp(['  cumulative iter: ', num2str(usr_par.cumulative_iter)]);
    disp '=============================';
    
    
    % change freq-dependend info @ usr_par.
    usr_par.whichFrq = ifreq;
    usr_par.sEventInfo = sObsPerFreq(ifreq).sEventInfo;
    usr_par.sEventObs  = sObsPerFreq(ifreq).sEventObs;
    
    if usr_par.cumulative_iter == 1
        stap = InvProps.stepInit;
    else
        stap = 1.0;
    end
    
    % prepare L-BFGS options
    options.verbose = true;
    options.max_iterations = cfe;
    options.init_step_length = stap;
    options.tolerance = 1e-13;
    options.output_file = output_log;
%     options.wolfe_try_to_increase_step_length = true;
%     % retrieve model
%     m = map_parameters_to_m(usr_par.Model(end), usr_par);
    
    % run L-BFGS
    [flag, mfinal, usr_par]=optlib_lbfgs(m, options, usr_par);
    
    m = mfinal;
end


