% travel time kernel calculation routine

% MAKE SURE V_OBS IS PRESENT!!!
[v_obs,~,~,~,~,~]=run_forward;

% run forward
cd ../code/
[v_rec,t,u_fw,v_fw,rec_x,rec_z]=run_forward;


% check what the seismograms look like
cd ../tools/
v_rec_3 = cat(3, [v_rec.x], [v_rec.y], [v_rec.z]);
plot_seismograms(v_rec_3,t,'velocity');

% make adjoint sources
adstf = make_all_adjoint_sources(v_rec,v_obs,t,'waveform_difference');

% check the adjoint source time functions
plot_vrec_to_adjointstf(t,v_rec.x,squeeze(adstf(1,:,:)));
plot_vrec_to_adjointstf(t,v_rec.y,squeeze(adstf(2,:,:)));
plot_vrec_to_adjointstf(t,v_rec.z,squeeze(adstf(3,:,:)));

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

% calculate the step length

% smooth the model