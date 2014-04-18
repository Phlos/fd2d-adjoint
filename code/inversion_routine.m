% travel time kernel calculation routine

% MAKE SURE V_OBS IS PRESENT!!!
cd ../code/
[v_obs,~,~,~,~,~]=run_forward;
save('../output/v_obs', 'v_obs', '-v7.3');

% run forward in test model
[v_iter(1),t,u_fw,v_fw,~,~]=run_forward;


% check what the seismograms look like
cd ../tools/
v_rec_3 = cat(3, [v_iter(1).x], [v_iter(1).y], [v_iter(1).z]);
plot_seismograms(v_rec_3,t,'velocity');
v_obs_3 = cat(3, [v_obs.x], [v_obs.y], [v_obs.z]);
plot_seismograms(v_obs_3,t,'velocity');

% make adjoint sources
cd ../tools
[adstf, misfit_iter(1)] = make_all_adjoint_sources(v_iter(1),v_obs,t,'waveform_difference','auto');
% adstf_iter(1).x = adstf(1,:,:);
% adstf_iter(1).y = adstf(2,:,:);
% adstf_iter(1).z = adstf(3,:,:);

% check the adjoint source time functions
plot_vrec_to_adjointstf(t,v_iter(1).x,squeeze(adstf(1,:,:)));
plot_vrec_to_adjointstf(t,v_iter(1).y,squeeze(adstf(2,:,:)));
plot_vrec_to_adjointstf(t,v_iter(1).z,squeeze(adstf(3,:,:)));


% run adjoint 
cd ../code/
K = run_adjoint(u_fw,v_fw,adstf,'waveform_difference');

% plot the kernels
cd ../tools/
[K, K_rel] = calculate_other_kernels(K);
plot_kernels_rho_mu_lambda_relative(K_rel);
print(gcf,'-dpng','-r400','../output/kernels-relative_rho-mu-lambda.png');
plot_kernels_rho_vs_vp_relative(K_rel);
print(gcf,'-dpng','-r400','../output/kernels-relative_rho-vs-vp.png');
plot_kernels_rho_mu_kappa_relative(K_rel);
print(gcf,'-dpng','-r400','../output/kernels-relative_rho-mu-kappa.png');

% calculate the step length and model update
[step,steplnArray,misfitArray] = calculate_step_length(1.5e21,2,misfitVrec,K,v_obs);
Params = update_model(K,step);
clearvars('u_fw');
clearvars('v_fw');
[v_iter(2),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Params);

% apply hard constraints
% -> no negative velocities
% -> 