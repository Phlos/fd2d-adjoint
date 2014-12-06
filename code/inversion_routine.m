
% preparation
path(path,'../input');
path(path,'../code');
path(path,'../code/propagation');
path(path,'../tools');
path(path,'../tools/misfits');
path(path,'../quivers');
path(path,'../mtit');

% number of iterations
niter = 2;
istart = 1;



% obtain useful parameters from input_parameters
[project_name, axrot, apply_hc, use_grav, parametrisation, ...
 rec_g, X, Z, normalise_misfits, stepInit] = get_input_info;




% MAKE SURE V_OBS IS PRESENT!!! AND SAVED!!!
 %[Model_real, v_obs, t_obs, props_obs, g_obs] = prepare_obs(modelnr);

% saving the observed variables
disp 'saving obs variables to matfile...'
savename = ['../output/',project_name,'.obs.all-vars.mat'];
save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp(['INVERSION RUN ', project_name]);
disp '-- using seis? ... YES!!'
if strcmp(use_grav , 'yes')
    disp '-- using grav? ... YES!!'
else
    disp '-- using grav? ... no'
end
if strcmp(apply_hc , 'yes')
    disp '-- using h.c.? ... YES!!'
else
    disp '-- using h.c.? ... no'
end
disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';


%=========================================================================



% save initial rho mu lambda (from input_parameters) as iter 1 Params values:
Model_start = update_model();

% % set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));

% plot very initial model
fig_mod = plot_model(Model_start,middle,parametrisation);
titel = [project_name,': model of iter 0'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter0.model.',parametrisation,'.png'];
print(fig_mod,'-dpng','-r400',figname);
close(fig_mod);

% apply hard constraints in initial model
if(strcmp(apply_hc,'yes'))
    % -> no negative velocities
    % -> mass of the Earth and/or its moment of inertia
    switch parametrisation
        case 'rhomulambda'
            [Model(1).rho, fig_hcupdate,~,~] = ...
                apply_hard_constraints(props_obs, Model_start.rho,axrot);
        case 'rhovsvp'
            Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model_start);
            [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
            Model(1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
        otherwise
            error('the parametrisation of the inversion was not recognised')
    end
    figname = ['../output/iter0.hard-constraints-rhoupdate.png'];
    titel = [project_name,': hc update of model 0'];
    mtit(fig_hcupdate, titel, 'xoff', 0.001, 'yoff', 0.05);
    print(fig_hcupdate,'-dpng','-r400',figname);
    close(fig_hcupdate);
    clearvars fig_rhoupdate
else
    %             disp 'model 1 is model start'
    Model(1) = Model_start;
end


for i = istart : niter;
%   if i > 1
        cd ../code;
        
        disp  ' ';
        disp  ' ';
        disp  ' ';
        disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
        disp '======================================';
        disp(['STARTING ITER ', num2str(i), ' OUT OF ', num2str(niter)]);
        disp '======================================';
        disp ' ';
        
        %% MODEL
        
        % plot model
        fig_mod = plot_model(Model(i),middle,parametrisation);
        titel = [project_name,': model of iter ', num2str(i)];
        mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = ['../output/iter',num2str(i),'.model.',parametrisation,'.png'];
        print(fig_mod,'-dpng','-r400',figname);
        close(fig_mod);
        clearvars('fig_mod');
        
        %% MISFITS
        
%         if strcmp(use_grav,'yes')
            % gravity field of current model
            [g(i), fig_grav] = calculate_gravity_field(Model(i).rho, rec_g);
%             figname = ['../output/iter',num2str(i),'.gravity_recordings.png'];
%             titel = [project_name,': gravity field of ', num2str(i), 'th model'];
%             mtit(fig_grav, titel, 'xoff', 0.001, 'yoff', 0.00001);
%             print(fig_grav, '-dpng', '-r400', figname);
            close(fig_grav);
            % comparison to real model:
            fig_grav_comp = plot_gravity_quivers(rec_g, g(i), g_obs, X, Z, Model(i).rho);
            figname = ['../output/iter',num2str(i),'.gravity_difference.png'];
            titel = [project_name,': gravity diff iter ', num2str(i), ' - real model'];
            mtit(fig_grav_comp, titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_grav_comp, '-dpng', '-r400', figname);
            close(fig_grav_comp);
            clearvars('fig_mod');
            
            %- calculate gravity misfit:
            % normalise the total misfit by first misfit (if already calc'd)
            if ( strcmp(normalise_misfits, 'byfirstmisfit')  && i==1)
                norm_misf_g = NaN;
            elseif (strcmp(normalise_misfits, 'byfirstmisfit') && i>1)
                norm_misf_g = misfit_g(1).total;
            else
                norm_misf_g = 1;
            end
            [g_src, misfit_g(i)] = make_gravity_sources(g(i), g_obs, norm_misf_g);
            clearvars norm_misf;
            
            % gravity misfit
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = misfit_g(i).total / sumgobs;
            disp ' ';
            disp(['gravity misfit for iter ',num2str(i),': ', ...
                num2str(misfit_g(i).total,'%3.2e')])
            disp(['   fraction of g_obs ',num2str(i),': ', ...
                num2str(div_by_gobs,'%3.2e')])
            if strcmp(normalise_misfits,'byfirstmisfit')
                disp(['   fraction of 1st misfit: ', ...
                    num2str(misfit_g(i).normd,'%3.2e')])
            end
            disp ' ';
%         else
%             misfit_g(i).total = NaN;
%             misfit_g(i).normd = NaN;
%         end

        

        
        %% SEISMIC MISFIT
        
        % run forward wave propagation
        disp ' ';
        disp(['iter ',num2str(i),': calculating forward wave propagation']);
%         clearvars u_fw v_fw;
        [v_iter(i),t,u_fw,v_fw,rec_x,rec_z]=run_forward(Model(i));
        close(clf);
        close(clf);
        close(clf);
        
        % make adjoint sources
        cd ../tools
        disp ' ';
        disp(['iter ',num2str(i),': calculating adjoint stf']);
        
        % normalise the total misfit by first misfit (if already calc'd)
        if ( strcmp(normalise_misfits, 'byfirstmisfit')  && i==1)
            norm_misf_s = NaN;
        elseif (strcmp(normalise_misfits, 'byfirstmisfit') && i>1)
            norm_misf_s = misfit_seis(1).total;
        else
            norm_misf_s = 1;
        end
        [adstf, misfit_seis(i)] = make_all_adjoint_sources( ...
                v_iter(i),v_obs,t,'waveform_difference','auto', norm_misf_s);
        clearvars norm_misf
        
       
        % plot seismogram difference
        fig_seisdif = plot_seismogram_difference(v_obs,v_iter(i),t);
        titel = [project_name,': difference between seismograms iter ', num2str(i), ' and obs'];
        mtit(fig_seisdif, titel, 'xoff', 0.001, 'yoff', 0.02);
        figname = ['../output/iter',num2str(i),'.seisdif.png'];
        print(fig_seisdif,'-dpng','-r400',figname);
        close(fig_seisdif);
        
%  end
        %% total misfit
        
        disp ' ';
        disp '=========================================='
        disp(['           MISFIT ITER ',num2str(i)]);
        
        if strcmp(use_grav,'yes')
            % gravity misfit
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = misfit_g(i).total / sumgobs;
            disp(['GRAVITY MISFIT FOR ITER ',num2str(i),': ', ...
                num2str(misfit_g(i).total,'%3.2e')])
            disp(['   fraction of g_obs:       ', ...
                num2str(div_by_gobs,'%3.2e')])
            if strcmp(normalise_misfits,'byfirstmisfit')
                disp(['   fraction of 1st misfit:  ', ...
                    num2str(misfit_g(i).normd,'%3.2e')])
            end
            disp ' ';
        end
        
        % seismic misfit
        sumvobs = sum(v_obs.x(:) .^2) + sum(v_obs.z(:) .^2);
        div_by_vobs = misfit_seis(i).total / sumvobs;
        disp ' ';
        disp(['SEISMIC MISFIT FOR ITER ',num2str(i),': ', ...
            num2str(misfit_seis(i).total,'%3.2e')])
        disp(['   fraction of v_obs:       ', ...
            num2str(div_by_vobs,'%3.2e')])
        if strcmp(normalise_misfits,'byfirstmisfit')
        disp(['   fraction of 1st misfit:  ', ...
            num2str(misfit_seis(i).normd,'%3.2e')])
        end
        disp ' ';
        
        if strcmp(use_grav,'yes')
            misfit(i) = misfit_seis(i).normd + misfit_g(i).normd;
        else
            misfit(i) = misfit_seis(i).normd;
        end
        
        disp(['TOTAL MISFIT FOR ITER ',num2str(i),': ', ...
            num2str(misfit(i),'%3.2e')])
        disp ' ';
        disp '=========================================='
        
        if i>1
        modeldifn(i) =   norm( (Model(i).rho(:) - Model(i-1).rho(:)) ./ Model(i-1).rho(:) ) ...
                          + norm( (Model(i).mu(:)  - Model(i-1).mu(:))  ./ Model(i-1).mu(:) ) ...
                          + norm( (Model(i).lambda(:) - Model(i-1).lambda(:)) ./ Model(i-1).lambda(:) );
        else
            modeldifn(i) = NaN;
        end 
        % plot misfit evolution
        fig_misfit = plot_misfit_evolution(misfit_seis,misfit_g,misfit,modeldifn);
        figname = '../output/misfit-evolution.png';
        mtit(fig_misfit, project_name, 'xoff', 0.001, 'yoff', 0.04);
        print(fig_misfit,'-dpng','-r400',figname);
        close(fig_misfit);
        clearvars fig_misfit
        
%  end
        %% GRAVITY KERNEL
        
        if strcmp(use_grav,'yes')
            
            % kernels only to be calculated when a next iteration will take place.
            if(i < niter)
                %- calculate gravity kernels
                disp ' ';
                disp(['iter ',num2str(i),': calculating gravity kernel']);
                % normalising the gravity kernel
                if strcmp(normalise_misfits,'byfirstmisfit')
                    normKg = 1.0 / misfit_g(1).total;
                else
                    normKg = 1.0;
                end
                % calculating the gravity kernel
                if i == 1
                    [Kg{i}, fig_Kg] = compute_kernels_gravity(g_src,rec_g,normKg,'no'); % 'no' is for plotting gravity kernel update
                else
                    [Kg{i}, fig_Kg] = compute_kernels_gravity(g_src,rec_g,normKg,'no'); % 'no' is for plotting gravity kernel update
                end
                
                %  plot gravity kernel
                figname = ['../output/iter',num2str(i),'.kernel_grav.rho.png'];
                titel = [project_name,': gravity kernel for iter ',num2str(i)];
                mtit(fig_Kg,titel, 'xoff', 0.001, 'yoff', 0.00001);
                print(fig_Kg,'-dpng','-r400',figname);
                close(fig_Kg);
                clearvars fig_Kg;
            end
        end
        
        
        %% SEISMIC KERNEL
        
        % kernels only to be calculated when a next iteration will take place.
        if(i < niter)
            %  if (i>1)
            % normalising the gravity kernel
            if strcmp(normalise_misfits,'byfirstmisfit')
                normKseis = 1.0 / misfit_seis(1).total;
            else
                normKseis = 1.0;
            end
            
            % run adjoint to obtain seismic kernels
            disp ' ';
            disp(['iter ',num2str(i),': calculating adjoint wave propagation']);
            cd ../code/
            Kseis(i) = run_adjoint(u_fw,v_fw,adstf,'waveform_difference',Model(i),normKseis);
            % storing kernels is not really necessary when allvars are saved at the end
%             disp 'storing kernels...'
%             kernelsavename = ['../output/iter', num2str(i),'.kernels.mat'];
%             save(kernelsavename,'Kseis', 'Kg','-v7.3');
            
            % empty the big variables so that the computer doesn't slow down.
            clearvars('u_fw');
            clearvars('v_fw');
%  end
%  disp 'hellooooo'
            % plot the kernels
            disp ' ';
            disp(['iter ',num2str(i),': plotting kernels']);
            cd ../tools/
            
            switch parametrisation
                case 'rhomulambda'
                    [~, K_reltemp] = calculate_other_kernels(Kseis(i), Model(i));
                    fig_knl = plot_kernels_rho_mu_lambda_relative(K_rel);
                    figname = ['../output/iter',num2str(i),'.kernels.relative.rho-mu-lambda.png'];
                    titel = [project_name,': seismic kernels (relative rhomulambda) for iter ',num2str(i)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    print(fig_knl,'-dpng','-r400',figname);
                    close(fig_knl);
                case 'rhovsvp'
                    [Kabs, K_reltemp] = calculate_other_kernels(Kseis(i), Model(i));
                    %                 K_reltemp = change_parametrisation_kernels('rhomulambda','rhovsvp',K_rel, Model(i));
                    fig_knl = plot_kernels_rho_vs_vp_relative(K_reltemp);
                    titel = [project_name,': seismic kernels (relative rhovsvp) for iter ',num2str(i)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(i),'.kernels.relative.rho-vs-vp.png'];
                    print(fig_knl,'-dpng','-r400',figname);
                    close(fig_knl);
                    fig_knl = plot_kernels_rho_vs_vp(Kabs);
                    titel = [project_name,': seismic kernels (absolute rhovsvp) for iter ',num2str(i)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(i),'.kernels.abs.rho-vs-vp.png'];
                    print(fig_knl,'-dpng','-r400',figname);
                    close(fig_knl);
                otherwise
                    error('unrecognised parametrisation for kernel plot');
            end
        end
        clearvars K_reltemp fig_knl;
        
%     disp 'hellooooo'
%     end
    
    
    %% COMBINE KERNELS & UPDATE MODEL
    
    
    % only update the model if we're going to a next model
    if (i<niter)
%    if i>1

        if strcmp(use_grav,'yes')
            % determine weight of relative kernels
            w_Kseis = 1;
            w_Kg = 1;
            
            disp(['seismic kernel 98th prctile:        ',num2str(prctile(abs(Kseis(i).rho.total(:)),98))]);
            disp(['gravity kernel 98th prctile:        ',num2str(prctile(abs(Kg{i}(:)),98))]);
            
            
            verhouding(i) = prctile(abs(Kseis(i).rho.total(:)),98) / prctile(abs(Kg{i}(:)),98);
            disp(['the ratio of seis and grav kernels: ',num2str(verhouding(i))]);
            disp(['the ratio of grav and seis weights: ',num2str(w_Kg/w_Kseis)]);
            
            % combine seismic and gravity kernels
            disp ' ';
            disp(['iter ',num2str(i),': combining gravity and seismic kernels']);
            Ktest = change_parametrisation_kernels('rhomulambda',parametrisation,Kseis(i),Model(i));
            Ktest.rho2.total = w_Kseis * Ktest.rho2.total  +  w_Kg * Kg{i};
            Ktest1 = change_parametrisation_kernels(parametrisation,'rhomulambda', Ktest,Model(i));
            K_total(i).rho.total    = Ktest1.rho.total;
            K_total(i).mu.total     = Ktest1.mu.total;
            K_total(i).lambda.total = Ktest1.lambda.total;
            clearvars('Ktest', 'Ktest1');
        else
            K_total(i) = Kseis(i);
        end
%     end

        % calculate the step length and model update
        disp ' ';
        disp(['iter ',num2str(i),': calculating step length']);
        if i==1;
            % basis for new step is initial step length
            [step(i), fig_lnsrch] = calculate_step_length(stepInit,i, ...
                misfit(i), misfit_seis, misfit_g, ...
                Model(i),K_total(i),v_obs, g_obs);
        elseif i>1;
            % basis for new step length is previous one.
            [step(i), fig_lnsrch] = calculate_step_length(step(i-1),i, ...
                misfit(i), misfit_seis, misfit_g, ...
                Model(i),K_total(i),v_obs, g_obs);
        end        
        % save linesearch figure
        figname = ['../output/iter',num2str(i),'.step-linesearch.png'];
        titel = [project_name,': linesearch for iter ',num2str(i)];
        mtit(fig_lnsrch,titel, 'xoff', 0.001, 'yoff', 0.00001);
        print(fig_lnsrch,'-dpng','-r400',figname);
        close(fig_lnsrch);
        clearvars fig_lnsrc;

        % actual model update
        disp ' ';
        disp(['iter ',num2str(i),': updating model']);
        Model(i+1) = update_model(K_total(i),step(i),Model(i),parametrisation);

        
        
        %% HARD CONSTRAINTS
        % apply hard constraints
        if(strcmp(apply_hc,'yes'))
            % -> no negative velocities
            % -> mass of the Earth and/or its moment of inertia
            switch parametrisation
                case 'rhomulambda'
                    [Model(i+1).rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model(i+1).rho,axrot);
                case 'rhovsvp'
                    Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model(i+1));
                    [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                    Model(i+1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
                otherwise
                    error('the parametrisation of the inversion was not recognised')
            end
            figname = ['../output/iter',num2str(i),'.hard-constraints-rhoupdate.png'];
            titel = [project_name,': hard constraints update for iter ',num2str(i)];
            mtit(fig_hcupdate,titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_hcupdate,'-dpng','-r400',figname);
            close(fig_hcupdate);
            clearvars fig_rhoupdate
        end
        
%         %% calculating model norm
%         modeldifnorm(i) =   norm( (Model(i+1).rho(:) - Model(i).rho(:)) ./ Model(i).rho(:) ) ...
%                           + norm( (Model(i+1).mu(:)  - Model(i).mu(:))  ./ Model(i).mu(:) ) ...
%                           + norm( (Model(i+1).lambda(:) - Model(i).lambda(:)) ./ Model(i).lambda(:) );
        
    end
    
    %% OUTPUT:
    
    % saving current variables to file (crash safeguard)
    disp 'saving all current variables...'
    clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
        'filenm_old', 'filenm_new', 'fig_knl');
%     exclude_vars = {'u_fw', 'v_fw'};
    savename = ['../output/',project_name,'.all-vars.mat'];
    save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');
    
end







%% WRAP-UP: misfit evolution & saving all variables to file

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|         ...FINISHING UP...         |';
disp '======================================';

% disp 'saving misfit evolution...'
% savename = ['../output/',project_name,'.misfit-evolution.mat'];
% save(savename, 'misfit', '-v7.3');

disp 'saving all current variables...'
clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
    'filenm_old', 'filenm_new', 'fig_knl');
savename = ['../output/',project_name,'.all-vars.mat'];
save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|               DONE!                |'
disp '======================================';
disp ' ';
disp(['  (this was inversion ',project_name,')']);