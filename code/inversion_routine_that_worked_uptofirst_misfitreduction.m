% inversion routine that worked up to the 1st real model update + misfit
% reduction!

load ../output/rho_anomaly_v_obs % I still had a v_obs for the rho anomaly in model 13

% run forward calculation with test model 11 (homogeneous)
[v_iter(1),t,u_fw,v_fw,~,~]=run_forward;

% look at the seismograms for a bit (not really vital)
cd ../tools/
v_rec_3 = cat(3, [v_iter(1).x], [v_iter(1).y], [v_iter(1).z]);
plot_seismograms(v_rec_3,t,'velocity');
v_obs_3 = cat(3, [v_obs.x], [v_obs.y], [v_obs.z]);
plot_seismograms(v_obs_3,t,'velocity');

% make adjoint sources
[adstf, misfit] = make_all_adjoint_sources(v_iter(1),v_obs,t,'waveform_difference','auto');

% check the adjoint source time functions (not really vital either)
plot_vrec_to_adjointstf(t,v_iter(1).x,squeeze(adstf(1,:,:)));

% run adjoint
cd ../code/
K = run_adjoint(u_fw,v_fw,adstf,'waveform_difference');

% produce the kernels & view them
cd ../tools/
[K, K_rel] = calculate_other_kernels(K);
plot_kernels_rho_mu_lambda_relative(K_rel);

% save the adstf and misfit from the previous run so they don't get
% overwritten
adstf_iter1 = adstf;
misfit_iter(1) = misfit;

% calculate the step length (this takes LONG).
[step,steplnArray,misfitArray] = calculate_step_length(1.5e21,1,misfit_iter(1),K,v_obs);

% update rho mu lambda model parameters
Params = update_model(K,step);

% clear big vars so that the new forward and adjoint calcs are not hindered
% by it (could've done this earlier)
clearvars('u_fw');
clearvars('v_fw');

% run the forward calculation with the updated model parameters
[v_iter(2),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Params);

% calculate the misfit of the 2nd inversion iteration
[adstf, misfit_iter(2)] = make_all_adjoint_sources(v_iter(2),v_obs,t,'waveform_difference','auto');

% the result: the misfit has decreased! Oh BOY.