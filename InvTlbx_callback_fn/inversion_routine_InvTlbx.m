% inversion routine that interacts with Christian's Optimisation Toolbox 
% code

%% preparation

% obtain useful parameters from input_parameters
[project_name, ~, apply_hc, use_grav, use_seis, fix_velocities, fix_density, ...
    use_matfile_startingmodel, starting_model, bg_model_type,...
    true_model_type, ~, change_freq_every, ...
    parametrisation, param_plot, rec_g, ~, ~, ~, ...
    normalise_misfits, InvProps.stepInit, smoothgwid] = get_input_info;

% kernels to be added in same parametrisation as inversion.
param_addknls = parametrisation;

% project folder
output_path = ['./output/',project_name,'/'];
mkdir('./output/',project_name)

% save input_parameters to this path
savename = [output_path,project_name,'.input_parameters.m'];
copyfile('./input/input_parameters.m',savename)

%% welcome
% 
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp(['INVERSION RUN ', project_name]);
disp(['-- parametrisation:  ', parametrisation]);
if strcmp(use_seis , 'yesseis')
    disp '-- using seis? ..... YES!!'
else
    disp '-- using seis? ..... no'
end
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
if strcmp(fix_density , 'yes')
    disp '-- fixing density? .... YES!!'
else
    disp '-- fixing density? .... no'
end
disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';


%% OBS

obs_file = [output_path,'/obs.all-vars.mat']
% freq consists of freqs.v_obs and freqs.frequency
if ((~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var'))  && exist(obs_file, 'file'))
    disp 'loading OBS file...'
    load(obs_file);
elseif (~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var'))
%    if istart == 1
        disp 'no OBS present, preparing obs...';
        [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path,true_model_type);
        save(obs_file, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');
%    else
%        error('iter > 1 but there are no observed properties!')
%    end
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

% % make initial model (from input_parameters)
% Model_start = update_model();

% make background model (for plotting purposes)
Model_bg = update_model(bg_model_type);

% % set the background values for plot_model to mode of the initial model
% middle.rho    = mode(Model_start.rho(:));
% middle.mu     = mode(Model_start.mu(:));
% middle.lambda = mode(Model_start.lambda(:));


% 1st model in a matfile?
% if istart == 1
    
    % if 1st iter model @ matfile, load matfile
    if strcmp(use_matfile_startingmodel,'yes')
        startmod = load(starting_model);
        Model_start = startmod.Model;
        fig_mod_start = plot_model_diff(Model_start, Model_bg, 'rhovsvp');
        close(fig_mod_start);
        clearvars startmod fig_mod_start;
    else
        Model_start = update_model();
    end
  
    
% end

%% calculate initial model misfit per frequency (if not present yet)

init_misfit_file = [output_path,'/initial_misfits.mat'];
% skip calculating initial misfits if they are already in the workspace
% or if there's a initial misfit file in the inversion directory
if ~strcmp(normalise_misfits, 'no')
    if  (~exist('misfit_init', 'var') && exist(init_misfit_file, 'file'))
        disp 'loading initial misfit file'
        load(init_misfit_file);
    elseif ~exist('misfit_init', 'var')
        disp 'calculating initial misfits'
        [misfit_init, misfit_init_traces] = calc_initial_misfits(Model_start, sObsPerFreq, g_obs);
        save(init_misfit_file, 'misfit_init', 'misfit_init_traces',  '-v6');
    else
        disp 'initial misfits already present... proceeding...';
    end
else
    misfit_init.total = NaN;
    misfit_init.seis = NaN; 
    misfit_init.grav = NaN;
end

%% set initial inversion values

% set the stage for multiple frequency inv.
whichFrq = 1;
cfe = change_freq_every;
cumulative_iter = 0;
nfreq = numel(sObsPerFreq);

if strcmp(use_seis, 'yesseis')
    sEventInfo = sObsPerFreq(whichFrq).sEventInfo;
    sEventObs   = sObsPerFreq(whichFrq).sEventObs;
else
    sEventInfo = 'noSeismicInfoUsed';
    sEventObs = 'noSeismicInfoUsed';
end




%% preparing input for the Inversion Toolbox

% output log
output_log = [output_path,'/lbfgs_output_log.txt'];

% prepare usr_par
usr_par.InvProps        = InvProps;
usr_par.Model_start     = Model_start;
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
usr_par.cumulative_iter = cumulative_iter;
if exist('Model_real', 'var')
    usr_par.Model_real  = Model_real;
end

% starting model
m = map_parameters_to_m(Model_start, usr_par);

% remove old fwd fields
TempFolder = [output_path,'/fwd_temp'];
if (exist(TempFolder, 'dir'))
    rmdir(TempFolder, 's');
end


%% calling onto the Inversion Toolbox

% loop over L-BFGS calls at different frequencies
for ifreq = 1:nfreq
   
    disp '=============================';
    disp '  going to a new frequency   ';
    disp(['  (frequency ',num2str(ifreq),'/',num2str(nfreq)',')']);
    disp(['  cumulative iter: ', num2str(usr_par.cumulative_iter + 1)]);
    disp '=============================';
    
    
    % change freq-dependend info @ usr_par.
    if strcmp(use_seis, 'yesseis')
        usr_par.whichFrq = ifreq;
        usr_par.sEventInfo = sObsPerFreq(ifreq).sEventInfo;
        usr_par.sEventObs  = sObsPerFreq(ifreq).sEventObs;
    end
    
    
    if usr_par.cumulative_iter == 0
        stap = InvProps.stepInit;
    else
        stap = 1.0;
    end
    
    % prepare L-BFGS options
    options.verbose = true;
    options.max_iterations = cfe;
    options.init_step_length = stap;
    options.grad_step_length = InvProps.stepInit;
    options.tolerance = 1e-13;
    options.output_file = output_log;
%     options.wolfe_try_to_increase_step_length = true;
%     % retrieve model
%     m = map_parameters_to_m(usr_par.Model(end), usr_par);
    
    % run L-BFGS (this is the actual inversion of all the iterations for
    % this frequency band)
    [flag, mfinal, usr_par]=optlib_lbfgs(m, options, usr_par);

    % set cumulative iter to itself - 1 because we're going to recalculate
    % the gradient of the last model at the new frequency
    usr_par.cumulative_iter = usr_par.cumulative_iter - 1;
    m = mfinal;
    
    % save optlib output for later use (if crash between freqs)
    save([output_path, '/freq-', num2str(ifreq, '%03d'), '.optlib_output.mat'], ...
        'ifreq', 'flag', 'm' , 'usr_par')
end

disp(['Finished inversion ', project_name, '!']);
