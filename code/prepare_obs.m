function [Model_real, v_obs, t_obs, props_obs, g_obs] = prepare_obs(varargin)

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

% real - starting model
Mstart = update_model(model_type);
fig_mod_diff = plot_model_diff(Model_real, Mstart, 'rhovsvp');
mtit(fig_mod_diff, 'real - starting model -- rho-vs-vp parametrisation');
figname = ['../output/obs.model-real-diff-starting.rhovsvp.png'];
print(fig_mod_diff, '-dpng', '-r400', figname);

% h.c. model properties for real model
props_obs = calculate_model_properties(Model_real.rho, 'x');

% gravity field of real model
[g_obs, fig_grav_obs] = calculate_gravity_field(Model_real.rho, rec_g);
figname = ['../output/obs.gravityrecordings.png'];
mtit(fig_grav_obs, 'gravity field of real model');
print(fig_grav_obs, '-dpng', '-r400', figname);

% real wave propagation
[v_obs,t_obs,~,~,~,~] = run_forward(Model_real);

v_0 = make_seismogram_zeros(v_obs);
fig_seis = plot_seismogram_difference(v_obs,v_0,t_obs,'nodiff');
titel = [project_name,': observed seismograms'];
mtit(fig_seis, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/obs.seis.png'];
print(fig_seis,'-dpng','-r400',figname);

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = ['../output/obs.all-vars.mat'];
save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

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