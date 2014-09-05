
% preparation
path(path,'../tools');
path(path,'../input');
path(path,'./propagation');

% number of iterations
niter = 7;

% initial step length;
stepInit = 3.5e14;

% obtain project name
 [project_name, axrot, apply_hc, parametrisation] = get_input_info;











% MAKE SURE V_OBS IS PRESENT!!! AND SAVED!!!
% [v_obs, t_obs, props_obs] = prepare_obs(modelnr);




% is v_obs saved yet?






% save initial rho mu lambda (from input_parameters) as iter 1 Params values:
Model(1) = update_model();


% % fig_mod = plot_model(Model(1));
% % figname = ['../output/iter',num2str(1),'.model.png'];
% % print(fig_mod,'-dpng','-r400',figname);
    
% % run forward in test model -- 1st iteration
% % disp 'calculating iter 1 forward'
% % [v_iter(1),t,u_fw,v_fw,~,~]=run_forward;
% 
% % % check what the seismograms look like
% % cd ../tools/
% % v_rec_3 = cat(3, [v_iter(1).x], [v_iter(1).y], [v_iter(1).z]);
% % plot_seismograms(v_rec_3,t,'velocity');
% % v_obs_3 = cat(3, [v_obs.x], [v_obs.y], [v_obs.z]);
% % plot_seismograms(v_obs_3,t,'velocity');

for i = 3:niter;
%  if i > 1
    cd ../code;
    
    disp  ' ';
    disp  ' ';
    disp  ' ';
    disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
    disp '======================================';
    disp(['STARTING ITER ', num2str(i)]);
    disp '======================================';
    disp ' ';
    
    % plot model


    switch parametrisation
        case 'rhomulambda'
            fig_mod = plot_model(Model(i));
            figname = ['../output/iter',num2str(i),'.model.rhomulambda.png'];
            print(fig_mod,'-dpng','-r400',figname);
            close(fig_mod);
        case 'rhovsvp'
            fig_mod = plot_model(Model(i),'rhovsvp');
            figname = ['../output/iter',num2str(i),'.model.rhovsvp.png'];
            print(fig_mod,'-dpng','-r400',figname);
            close(fig_mod);
        otherwise
            error('unrecognised parametrisation for model plot');
    end
    
    % run forward wave propagation 
    disp ' ';
    disp(['iter ',num2str(i),': calculating forward wave propagation']);
    [v_iter(i),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Model(i));
    close(clf);
    close(clf);
    close(clf);
    

    % make adjoint sources
    cd ../tools
    disp ' ';
    disp(['iter ',num2str(i),': calculating adjoint stf']);
    [adstf, misfit_iter(i)] = make_all_adjoint_sources(v_iter(i),v_obs,t,'waveform_difference','auto');
    
    % plot seismogram difference
    fig_seisdif = plot_seismogram_difference(v_obs,v_iter(i),t);
    figname = ['../output/iter',num2str(i),'.seisdif.png'];
    print(fig_seisdif,'-dpng','-r400',figname);
    close(fig_seisdif);
    
    disp ' ';
    disp(['MISFIT FOR ITER ',num2str(i),': ', ...
          num2str(misfit_iter(i).total,'%3.2e')])
    disp ' ';
    
    % run adjoint
    disp ' ';
    disp(['iter ',num2str(i),': calculating adjoint wave propagation']);
    cd ../code/
    K(i) = run_adjoint(u_fw,v_fw,adstf,'waveform_difference',Model(i));
    
    % empty the big variables so that the computer doesn't slow down.
    clearvars('u_fw');
    clearvars('v_fw');

    
    % plot the kernels
    disp ' ';
    disp(['iter ',num2str(i),': plotting kernels']);
    cd ../tools/
    [K_abs(i), K_rel(i)] = calculate_other_kernels(K(i), Model(i));
    switch parametrisation
        case 'rhomulambda'
            fig_knl = plot_kernels_rho_mu_lambda_relative(K_rel(i));
            figname = ['../output/iter',num2str(i),'.kernels.relative.rho-mu-lambda.png'];
            print(fig_knl,'-dpng','-r400',figname);
    close(fig_knl);
        case 'rhovsvp'
            fig_knl = plot_kernels_rho_vs_vp_relative(K_rel(i));
            figname = ['../output/iter',num2str(i),'.kernels.relative.rho-vs-vp.png'];
            print(fig_knl,'-dpng','-r400',figname);
            close(fig_knl);
        otherwise
            error('unrecognised parametrisation for kernel plot');
    end
%     
%     disp 'kernels rho mu kappa relative'
%     plot_kernels_rho_mu_kappa_relative(K_rel(i));
% %     figname = ['../output/kernels-relative_rho-mu-kappa_iter-',num2str(i),'.png'];
%     print(gcf,'-dpng','-r400',figname);
%     close(gcf);

%  end    
    % calculate the step length and model update
    disp ' ';
    disp(['iter ',num2str(i),': calculating step length']);
    if i==1;
        [step(i)] = calculate_step_length(stepInit,i,misfit_iter(i), ...
                                Model(i),K_rel(i),v_obs);
    elseif i>1;
        % basis for new step length is previous one.
        [step(i)] = calculate_step_length(step(i-1),i,misfit_iter(i), Model(i),K_rel(i),v_obs);
    end

    
    disp ' ';
    disp(['iter ',num2str(i),': updating model']);
    %     if i>1
    Model(i+1) = update_model(K_rel(i),step(i),Model(i));
    %     end
    
% end

    % apply hard constraints
    if(strcmp(apply_hc,'yes'))
        % -> no negative velocities
        % -> mass of the Earth and/or its moment of inertia
        switch parametrisation
            case 'rhomulambda'
                [Model(i+1).rho, fig_rhoupdate,~,~] = ...
                    apply_hard_constraints(props_obs, Model(i+1).rho,axrot);
            case 'rhovsvp'
                Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model(i+1));
                [Model_rhovsvp.rho, fig_rhoupdate,~,~] = ...
                    apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                Model(i+1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
            otherwise
                error('the parametrisation of the inversion was not recognised')
        end
        figname = ['../output/iter',num2str(i),'.hard-constraints-rhoupdate.png'];
        print(fig_rhoupdate,'-dpng','-r400',figname);
        close(fig_rhoupdate);
    end
    
    % OUTPUT:
    
    % save kernels per iter
    filenm_old = ['../output/', project_name, '.kernels.mat'];
    filenm_new = ['../output/iter', num2str(i),'.kernels.mat'];
    movefile(filenm_old, filenm_new);

    % save v_rec per iter
    filenm_old = ['../output/', project_name, '.v_rec.mat'];
    filenm_new = ['../output/iter', num2str(i),'.v_rec.mat'];
    movefile(filenm_old, filenm_new);
    
    % rename linesearch step figure
    filenm_old = '../output/iter.step-linesearch.png';
    filenm_new = ['../output/iter', num2str(i), '.step-linesearch.png'];
    movefile(filenm_old, filenm_new);
    
end

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp 'FINISHING UP WITH THE LAST MODEL...'
disp '======================================';

% PLOT MODEL
fig_mod = plot_model(Model(niter+1));
figname = ['../output/iter',num2str(niter+1),'.model.png'];
print(fig_mod,'-dpng','-r400',figname);
close(fig_mod);


% FORWARD PROPAGATION
disp(['iter ',num2str(niter+1),': calculating forward wave propagation']);
[v_iter(niter+1),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Model(niter+1));
close(clf);
close(clf);
close(clf);
% empty the big variables so that the computer doesn't slow down.
clearvars('u_fw');
clearvars('v_fw');

% save v_rec per iter
filenm_old = ['../output/', project_name, '.v_rec.mat'];
filenm_new = ['../output/iter', num2str(niter+1),'.v_rec.mat'];
movefile(filenm_old, filenm_new);


% MISFIT:
cd ../tools
disp(['iter ',num2str(niter+1),': calculating adjoint stf']);
[adstf, misfit_iter(niter+1)] = make_all_adjoint_sources(v_iter(niter+1),v_obs,t,'waveform_difference','auto');

% plot seismogram difference
fig_seisdif = plot_seismogram_difference(v_obs,v_iter(niter+1),t);
%     figname = ['../output/iter',num2str(i),'.seisdif-', num2str(misfit_iter(i).total, '%3.2e'), '.png'];
figname = ['../output/iter',num2str(niter+1),'.seisdif.png'];
print(fig_seisdif,'-dpng','-r400',figname);
close(fig_seisdif);

disp ' ';
disp(['MISFIT FOR ITER ',num2str(niter+1),': ', ...
      num2str(misfit_iter(niter+1).total,'%3.2e')])
disp ' ';

% plot misfit evolution
for i=1:niter+1
    misfit(i) = misfit_iter(i).total;
end

fig_misfit = figure;
subplot(1,2,1);
plot(misfit);
title('misfit evolution');
subplot(1,2,2);
semilogy(misfit);
title('misfit evolution - semilogarithmic in misfit');
% subplot(1,3,3);
% loglog(misfit);
figname = ['../output/misfit-evolution.png'];
print(fig_misfit,'-dpng','-r400',figname);
close(fig_misfit);

disp 'saving misfit evolution...'
savename = ['../output/',project_name,'.misfit-evolution.mat'];
save(savename, 'misfit', '-v7.3');

disp 'saving all current variables...'
clearvars('fig_misfit', 'figname', 'fig_seisdif', 'fig_mod', ...
          'filenm_old', 'filenm_new', 'fig_knl');
savename = ['../output/',project_name,'.all-vars.mat'];
save(savename);

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|               DONE!                |'
disp '======================================';