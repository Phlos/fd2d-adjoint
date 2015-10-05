% resume iterations at frequency startfreq

startfreq = 4;

disp ' ';
disp(['STARTFREQ = ', num2str(startfreq)]);
disp ' ';


%% initialise 

% initial stuff
% obtain useful parameters from input_parameters
[project_name, ~, apply_hc, use_grav, use_seis, fix_velocities, ...
    use_matfile_startingmodel, starting_model, bg_model_type,...
    true_model_type, ~, change_freq_every, ...
    parametrisation, param_plot, rec_g, ~, ~, ~, ...
    ~, InvProps.stepInit, smoothgwid] = get_input_info;
output_path = ['./output/',project_name,'/'];

% load stuff
disp 'loading stuff'
obs_file = [output_path,'/obs.all-vars.mat'];
load(obs_file);
init_misfit_file = [output_path,'/initial_misfits.mat'];
load(init_misfit_file);
load([output_path, '/freq-', num2str(startfreq-1, '%03d'), '.optlib_output.mat']);

cfe = change_freq_every;
nfreq = length(sObsPerFreq);
output_log = [output_path,'/lbfgs_output_log.txt'];

%% iterations
% start iterations all over
for ifreq = startfreq:nfreq
   
    disp '=============================';
    disp '  going to a new frequency   ';
    disp(['  (frequency ',num2str(ifreq),'/',num2str(nfreq)',')']);
    disp(['  cumulative iter: ', num2str(usr_par.cumulative_iter)]);
    disp '=============================';
    
    
    % change freq-dependend info @ usr_par.
    if strcmp(use_seis, 'yesseis')
        usr_par.whichFrq = ifreq;
        usr_par.sEventInfo = sObsPerFreq(ifreq).sEventInfo;
        usr_par.sEventObs  = sObsPerFreq(ifreq).sEventObs;
    end
    
    if usr_par.cumulative_iter == 1
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
    
    % run L-BFGS
    [flag, mfinal, usr_par]=optlib_lbfgs(m, options, usr_par);

    m = mfinal;
    
        % save optlib output for later use (if crash between freqs)
    save([output_path, '/freq-', num2str(ifreq, '%03d'), '.optlib_output.mat'], ...
        'ifreq', 'flag', 'm' , 'usr_par')
end