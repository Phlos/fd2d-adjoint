% travel time kernel calculation routine

% number of iterations
niter = 5;

% initial step length;
stepInit = 3.0e21;

% MAKE SURE V_OBS IS PRESENT!!!
% cd ../code/
% [v_obs,~,~,~,~,~]=run_forward;
% save('../output/v_obs', 'v_obs', '-v7.3');


% save initial rho mu lambda as iter 1 Params values:
% Params(1) = update_model();
    
% run forward in test model -- 1st iteration
% disp 'calculating iter 1 forward'
% [v_iter(1),t,u_fw,v_fw,~,~]=run_forward;

% % check what the seismograms look like
% cd ../tools/
% v_rec_3 = cat(3, [v_iter(1).x], [v_iter(1).y], [v_iter(1).z]);
% plot_seismograms(v_rec_3,t,'velocity');
% v_obs_3 = cat(3, [v_obs.x], [v_obs.y], [v_obs.z]);
% plot_seismograms(v_obs_3,t,'velocity');

for i = 1:niter;
    disp(['starting iter ', num2str(i)]);
    
    % run forward in subsequent iterations >1
    if i>1
        disp(['calculating iter ',num2str(i),' forward']);
        [v_iter(i),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Params(i));
    end
    
    % make adjoint sources
    cd ../tools
    disp(['calculating iter ',num2str(i),' stf']);
    [adstf, misfit_iter(i)] = make_all_adjoint_sources(v_iter(i),v_obs,t,'waveform_difference','auto');
    
    misfit_iter(i)
    % adstf_iter(1).x = adstf(1,:,:);
    % adstf_iter(1).y = adstf(2,:,:);
    % adstf_iter(1).z = adstf(3,:,:);
    

%     % check the adjoint source time functions
%     plot_vrec_to_adjointstf(t,v_iter(i).x,squeeze(adstf(1,:,:)));
%     plot_vrec_to_adjointstf(t,v_iter(i).y,squeeze(adstf(2,:,:)));
%     plot_vrec_to_adjointstf(t,v_iter(i).z,squeeze(adstf(3,:,:)));
    
    
    % run adjoint
    cd ../code/
    disp(['calculating iter ',num2str(i),' adjoint']);
    K(i) = run_adjoint(u_fw,v_fw,adstf,'waveform_difference',Params(i));
%     end % this one is only there because I already calculated the first adjoint
    
    % plot the kernels
    cd ../tools/
    [K_abs(i), K_rel(i)] = calculate_other_kernels(K(i));
    
    disp 'kernels rho mu lambda relative'
    plot_kernels_rho_mu_lambda_relative(K_rel);
    figname = ['../output/kernels-relative_rho-mu-lambda_iter-',num2str(i),'.png'];
    print(gcf,'-dpng','-r400',figname);
    close(gcf);
    
    disp 'kernels rho vs vp relative'
    plot_kernels_rho_vs_vp_relative(K_rel);
    figname = ['../output/kernels-relative_rho-vs-vp_iter-',num2str(i),'.png'];
    print(gcf,'-dpng','-r400',figname);
    close(gcf);
    
    disp 'kernels rho mu kappa relative'
    plot_kernels_rho_mu_kappa_relative(K_rel);
    figname = ['../output/kernels-relative_rho-mu-kappa_iter-',num2str(i),'.png'];
    print(gcf,'-dpng','-r400',figname);
    close(gcf);
%end % this one is only there because I am debugging in a stupid way - refers back to if i>2
    
    % calculate the step length and model update
    disp 'calculating step length...'
    if i==1;
        [step(i),steplnArray,misfitArray] = ...
          calculate_step_length(stepInit,i,misfit_iter(i), ...
                                Params(i),K_rel(i),v_obs);
    elseif i>1;
        % basis for new step length is previous one.
        [step(i),steplnArray,misfitArray] = ...
          calculate_step_length(step(i-1),i,misfit_iter(i), ...
                                Params(i),K_rel(i),v_obs);
    end

    
    
    % apply hard constraints
    % -> no negative velocities
    % -> mass of the Earth and/or its moment of inertia

    
    disp 'updating model'
    %     if i>1
    Params(i+1) = update_model(K_rel(i),step(i),Params(i));
%     plot_model;
    %     end
    
    % empty the big variables so that the computer doesn't slow down.
    clearvars('u_fw');
    clearvars('v_fw');
    

    
end

