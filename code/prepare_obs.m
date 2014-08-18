function [v_obs, t_obs, props_obs] = prepare_obs(modelnr)

mod_real = update_model(modelnr);
props_obs = calculate_model_properties(mod_real.rho, 'x');
[v_obs,t_obs,~,~,~,~] = run_forward(mod_real);

close all;

end