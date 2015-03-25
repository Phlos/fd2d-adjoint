function [Model_real, freq, t_obs, props_obs, g_obs] = prepare_obs(varargin)

input_parameters;

Model_real = checkargs(varargin(:));

% Model_real = update_model(modelnr);

% % setting plot_model background values for all models based on Tromp 2005
% trompmodels = [10:39];
% if any(modelnr==trompmodels)
%     disp 'real model based on Tromp'
%     middle_real.rml = [2600 2.66e10    2.42e10];
%     middle_real.rvv = [2600 3198.55736 5797.87759];
% end

% plotting the real model
fig_mod_real = plot_model(Model_real, 'rhovsvp');
mtit(fig_mod_real, 'real model -- rho-vs-vp parametrisation');
figname = ['../output/obs.model.rhovsvp.png'];
print(fig_mod_real, '-dpng', '-r400', figname);
close(fig_mod_real);

% real - starting model
Mstart = update_model(model_type);
fig_mod_diff = plot_model_diff(Model_real, Mstart, 'rhovsvp');
mtit(fig_mod_diff, 'real - starting model -- rho-vs-vp parametrisation');
figname = ['../output/obs.model-real-diff-starting.rhovsvp.png'];
print(fig_mod_diff, '-dpng', '-r400', figname);
close(fig_mod_diff);

% h.c. model properties for real model
props_obs = calculate_model_properties(Model_real.rho, 'x');

% gravity field of real model
[g_obs, fig_grav_obs] = calculate_gravity_field(Model_real.rho, rec_g);
figname = ['../output/obs.gravityrecordings.png'];
mtit(fig_grav_obs, 'gravity field of real model');
print(fig_grav_obs, '-dpng', '-r400', figname);
close(fig_grav_obs);

%% source frequency dependent
% real wave propagation
for ii = 1:length(f_maxlist)
    
    % make source-time function per frequency
    freq(ii).frequency = f_maxlist(ii);
    freq(ii).stf = make_stf_wrapperscript(freq(ii).frequency);
    [freq(ii).v_obs,t_obs,~,~,~,~] = run_forward(Model_real, freq(ii).stf);

    % run wave propagation per frequency
    v_0 = make_seismogram_zeros(freq(ii).v_obs);
    fig_seis = plot_seismogram_difference(freq(ii).v_obs,v_0,t_obs,'nodiff');
    titel = [project_name,': observed seismograms f_max = ', num2str(freq(ii).frequency), ' Hz'];
    mtit(fig_seis, titel, 'xoff', 0.001, 'yoff', -0.05);
    figname = ['../output/obs.seis.fmax-',num2str(freq(ii).frequency,'%.2e'),'.png'];
    print(fig_seis,'-dpng','-r400',figname);
    close(fig_seis);
end

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = ['../output/obs.all-vars.mat'];
save(savename, 'freq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

close all;

end

function Model_real = checkargs(args)

nargs = length(args);

if nargs ~= 1
    error('wrong input to prepare_obs !')
else
    if isnumeric(args{1})
        modelnr = args{1};
        Model_real = update_model(modelnr);
    elseif isstruct(args{1})
        Model_real = args{1};
    end
    
end


end