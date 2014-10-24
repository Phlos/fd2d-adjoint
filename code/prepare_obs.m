function [Model_real, v_obs, t_obs, props_obs, g_obs] = prepare_obs(modelnr)

input_parameters;

Model_real = update_model(modelnr);

% plotting the real model
fig_mod_real = plot_model(Model_real, 'rhovsvp');
mtit(fig_mod_real, 'real model - rho-vs-vp parametrisation');
figname = ['../output/obs.model.rhovsvp.png'];
print(fig_mod_real, '-dpng', '-r400', figname);

% h.c. model properties for real model
props_obs = calculate_model_properties(Model_real.rho, 'x');

% gravity field of real model
[g_obs, fig_grav_obs] = calculate_gravity_field(Model_real.rho, rec_g);
figname = ['../output/obs.gravityrecordings.png'];
mtit(fig_grav_obs, 'gravity field of real model');
print(fig_grav_obs, '-dpng', '-r400', figname);

% real wave propagation
[v_obs,t_obs,~,~,~,~] = run_forward(Model_real);

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = ['../output/obs.all-vars.mat'];
save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

close all;

end