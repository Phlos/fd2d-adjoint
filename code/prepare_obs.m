function [v_obs, t_obs, props_obs] = prepare_obs(modelnr)

Model_real = update_model(modelnr);

% plotting the real model
fig_mod_real = plot_model(Model_real, 'rhovsvp');
mtit(fig_mod_real, 'real model - rho-vs-vp parametrisation');
figname = ['../output/model.real.rhovsvp.png'];
print(fig_mod_real, '-dpng', '-r400', figname);

props_obs = calculate_model_properties(Model_real.rho, 'x');
[v_obs,t_obs,~,~,~,~] = run_forward(Model_real);

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = ['../output/obs.all-vars.mat'];
save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', '-v7.3');

close all;

end