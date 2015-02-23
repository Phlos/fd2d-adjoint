
% preparation
path(path,'../input');
path(path,'../code');
path(path,'../code/propagation');
path(path,'../tools');
path(path,'../tools/misfits');
path(path,'../quivers');
path(path,'../mtit');

% number of iterations
InvProps.niter = 1;
InvProps.istart = 1;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, axrot, apply_hc, use_grav, parametrisation, ...
 rec_g, X, Z, misfit_type, normalise_misfits, InvProps.stepInit] = get_input_info;

%% just to get 'middle' right
% save initial rho mu lambda (from input_parameters) as iter 1 Params values:
Model_start = update_model();

% % set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));


%% OBS

% MAKE SURE V_OBS IS PRESENT!!! AND SAVED!!!
%[Model_real, v_obs, t_obs, props_obs, g_obs] = prepare_obs(modelnr);
 


% saving the observed variables
disp 'saving obs variables to matfile...'
savename = ['../output/',project_name,'.obs.all-vars.mat'];
save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

%% INVERSION

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
set(fig_mod,'Renderer','painters');
titel = [project_name,': model of iter 0'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter0.model.',parametrisation,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
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


for iter = InvProps.istart : InvProps.niter;
%   if i > 1
        cd ../code;
        
        disp  ' ';
        disp  ' ';
        disp  ' ';
        disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
        disp '======================================';
        disp(['STARTING ITER ', num2str(iter), ' OUT OF ', num2str(InvProps.niter)]);
        disp '======================================';
        disp ' ';
        
        %% MODEL
        
        % plot model
        fig_mod = plot_model(Model(iter),middle,parametrisation);
        titel = [project_name,': model of iter ', num2str(iter)];
        mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = ['../output/iter',num2str(iter),'.model.',parametrisation,'.png'];
        print(fig_mod,'-dpng','-r400',figname);
        close(fig_mod);
        clearvars('fig_mod');
        
        %% misfit
        
%         if strcmp(use_grav,'yes')
            % gravity field of current model
            [g(iter), fig_grav] = calculate_gravity_field(Model(iter).rho, rec_g);
%             figname = ['../output/iter',num2str(i),'.gravity_recordings.png'];
%             titel = [project_name,': gravity field of ', num2str(i), 'th model'];
%             mtit(fig_grav, titel, 'xoff', 0.001, 'yoff', 0.00001);
%             print(fig_grav, '-dpng', '-r400', figname);
            close(fig_grav);
            % comparison to real model:
            fig_grav_comp = plot_gravity_quivers(rec_g, g(iter), g_obs, X, Z, Model(iter).rho);
            figname = ['../output/iter',num2str(iter),'.gravity_difference.png'];
            titel = [project_name,': gravity diff iter ', num2str(iter), ' - real model'];
            mtit(fig_grav_comp, titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_grav_comp, '-dpng', '-r400', figname);
            close(fig_grav_comp);
            clearvars('fig_mod');
            
            %- calculate gravity misfit:
            % normalise the total misfit by first misfit (if already calc'd)
            if ( strcmp(normalise_misfits, 'byfirstmisfit')  && iter==1)
                norm_misf_g = NaN;
            elseif (strcmp(normalise_misfits, 'byfirstmisfit') && iter>1)
                norm_misf_g = InvProps.misfit_g(1).total;
            else
                norm_misf_g = 1;
            end
            [g_src, InvProps.misfit_g(iter)] = make_gravity_sources(g(iter), g_obs, norm_misf_g);
            clearvars norm_misf;
            
            % gravity misfit
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = InvProps.misfit_g(iter).total / sumgobs;
            disp ' ';
            disp(['gravity misfit for iter ',num2str(iter),': ', ...
                num2str(InvProps.misfit_g(iter).total,'%3.2e')])
            disp(['   fraction of g_obs ',num2str(iter),': ', ...
                num2str(div_by_gobs,'%3.2e')])
            if strcmp(normalise_misfits,'byfirstmisfit')
                disp(['   fraction of 1st misfit: ', ...
                    num2str(InvProps.misfit_g(iter).normd,'%3.2e')])
            end
            disp ' ';
%         else
%             misfit_g(i).total = NaN;
%             misfit_g(i).normd = NaN;
%         end

        

        
        %% SEISMIC misfit
        
        % run forward wave propagation
        disp ' ';
        disp(['iter ',num2str(iter),': calculating forward wave propagation']);
%         clearvars u_fw v_fw;
        [v_iter{iter},t,u_fw,v_fw,rec_x,rec_z]=run_forward(Model(iter));
        close(clf);
        close(clf);
        close(clf);
        
        % plot seismogram difference
        fig_seisdif = plot_seismogram_difference(v_obs,v_iter{iter},t);
        titel = [project_name,': difference between seismograms iter ', num2str(iter), ' and obs'];
        mtit(fig_seisdif, titel, 'xoff', 0.001, 'yoff', 0.02);
        figname = ['../output/iter',num2str(iter),'.seisdif.png'];
        print(fig_seisdif,'-dpng','-r400',figname);
        close(fig_seisdif);
        
        % make adjoint sources
        cd ../tools
        disp ' ';
        disp(['iter ',num2str(iter),': calculating adjoint stf']);
        

        
        % really make adjoint sources
        [adstf{iter}, InvProps.misfit_seis{iter}] = calc_misfitseis_adstf(misfit_type,t,v_iter{iter},v_obs);
%         [adstf, InvProps.misfit_seis(i)] = make_all_adjoint_sources( ...
%                 v_iter{i},v_obs,t,'waveform_difference','auto', norm_misf_s);
        % determine misfit normalisation factor by first misfit (if already calc'd)
        if ( strcmp(normalise_misfits, 'byfirstmisfit') )
             InvProps.misfit_seis{iter}.normd = InvProps.misfit_seis{iter}.total / ...
                                    InvProps.misfit_seis{1}.total;
%             norm_misf_s = NaN;
%         elseif (strcmp(normalise_misfits, 'byfirstmisfit') && iter>1)
%             norm_misf_s = InvProps.misfit_seis(1).total;
        else
            InvProps.misfit_seis{iter}.normd = InvProps.misfit_seis{iter}.total;
%             norm_misf_s = 1;
        end
        
        
%         % plot stf to adstf for one station
%         stationnr = 5;
%         fig_stftoadstf = plot_seismogram_difference(stf, disp, vel, adstf, t);
%         titel = [project_name,': stf to adstf at station', num2str(stationnr), ', iter ', num2str(i), ' and obs'];
%         mtit(fig_stftoadstf, titel, 'xoff', 0.001, 'yoff', 0.02);
%         figname = ['../output/iter',num2str(i),'.stf-to-adstf-station',num2str(stationnr),'.png'];
%         print(fig_stftoadstf,'-dpng','-r400',figname);
%         close(fig_stftoadstf);

        
%  end
        %% total misfit
        
        disp ' ';
        disp '=========================================='
        disp(['           misfit ITER ',num2str(iter)]);
        
        if strcmp(use_grav,'yes')
            % gravity misfit
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = InvProps.misfit_g(iter).total / sumgobs;
            disp(['GRAVITY misfit FOR ITER ',num2str(iter),': ', ...
                num2str(InvProps.misfit_g(iter).total,'%3.2e')])
            disp(['   fraction of g_obs:       ', ...
                num2str(div_by_gobs,'%3.2e')])
            if strcmp(normalise_misfits,'byfirstmisfit')
                disp(['   fraction of 1st misfit:  ', ...
                    num2str(InvProps.misfit_g(iter).normd,'%3.2e')])
            end
            disp ' ';
        end
        
        % seismic misfit
        sumvobs = 0;
        for ii = 1:length(v_obs)
            comp = fieldnames(v_obs{ii});
            for icomp = 1:length(comp)
            sumvobs = sumvobs + sum(v_obs{ii}.(comp{icomp}) .^2);
            end
        end
        div_by_vobs = InvProps.misfit_seis{iter}.total / sumvobs;
        disp ' ';
        disp(['SEISMIC misfit FOR ITER ',num2str(iter),': ', ...
            num2str(InvProps.misfit_seis{iter}.total,'%3.2e')])
        disp(['   fraction of v_obs:       ', ...
            num2str(div_by_vobs,'%3.2e')])
        if strcmp(normalise_misfits,'byfirstmisfit')
        disp(['   fraction of 1st misfit:  ', ...
            num2str(InvProps.misfit_seis{iter}.normd,'%3.2e')])
        end
        disp ' ';
        
        if strcmp(use_grav,'yes')
            InvProps.misfit(iter) = InvProps.misfit_seis{iter}.normd + InvProps.misfit_g(iter).normd;
        else
            InvProps.misfit(iter) = InvProps.misfit_seis{iter}.normd;
        end
        
        disp(['TOTAL misfit FOR ITER ',num2str(iter),': ', ...
            num2str(InvProps.misfit(iter),'%3.2e')])
        disp ' ';
        disp '=========================================='
        
        if iter>1
        InvProps.modeldifn(iter) =   norm( (Model(iter).rho(:) - Model(iter-1).rho(:)) ./ Model(1).rho(:) ) ...
                          + norm( (Model(iter).mu(:)  - Model(iter-1).mu(:))  ./ Model(1).mu(:) ) ...
                          + norm( (Model(iter).lambda(:) - Model(iter-1).lambda(:)) ./ Model(1).lambda(:) );
        else
            InvProps.modeldifn(iter) = NaN;
        end 
        % plot misfit evolution
        fig_misfit = plot_misfit_evolution(InvProps);
        figname = '../output/misfit-evolution.pdf';
        mtit(fig_misfit, project_name, 'xoff', 0.001, 'yoff', 0.04);
        print(fig_misfit,'-dpdf','-r400',figname);
        close(fig_misfit);
        clearvars fig_misfit
        
%  end
        %% GRAVITY KERNEL
        
        if strcmp(use_grav,'yes')
            
            % kernels only to be calculated when a next iteration will take place.
            if(iter < InvProps.niter)
                %- calculate gravity kernels
                disp ' ';
                disp(['iter ',num2str(iter),': calculating gravity kernel']);
                % normalising the gravity kernel
                if strcmp(normalise_misfits,'byfirstmisfit')
                    normKg = 1.0 / InvProps.misfit_g(1).total;
                else
                    normKg = 1.0;
                end
                % calculating the gravity kernel
                if iter == 1
                    [Kg{iter}, fig_Kg] = compute_kernels_gravity(g_src,rec_g,normKg,'no'); % 'no' is for plotting gravity kernel update
                else
                    [Kg{iter}, fig_Kg] = compute_kernels_gravity(g_src,rec_g,normKg,'no'); % 'no' is for plotting gravity kernel update
                end
                
                %  plot gravity kernel
                figname = ['../output/iter',num2str(iter),'.kernel_grav.rho.png'];
                titel = [project_name,': gravity kernel for iter ',num2str(iter)];
                mtit(fig_Kg,titel, 'xoff', 0.001, 'yoff', 0.00001);
                print(fig_Kg,'-dpng','-r400',figname);
                close(fig_Kg);
                clearvars fig_Kg;
            end
        end
        
        
        %% SEISMIC KERNEL
        
        % kernels only to be calculated when a next iteration will take place.
        if(iter < InvProps.niter)
            %  if (i>1)
            % normalising the gravity kernel
            if strcmp(normalise_misfits,'byfirstmisfit')
                normKseis = 1.0 / InvProps.misfit_seis{1}.total;
            else
                normKseis = 1.0;
            end
            
            % run adjoint to obtain seismic kernels
            disp ' ';
            disp(['iter ',num2str(iter),': calculating adjoint wave propagation']);
            cd ../code/
            Kseis(iter) = run_adjoint(u_fw,v_fw,adstf{iter},'waveform_difference',Model(iter),normKseis);
            % storing kernels is not really necessary when allvars are saved at the end
%             disp 'storing kernels...'
%             kernelsavename = ['../output/iter', num2str(i),'.kernels.mat'];
%             save(kernelsavename,'Kseis', 'Kg','-v7.3');
            

%  end
%  disp 'hellooooo'
            % plot the kernels
            disp ' ';
            disp(['iter ',num2str(iter),': plotting kernels']);
            cd ../tools/
            
            switch parametrisation
                case 'rhomulambda'
                    [~, K_reltemp] = calculate_other_kernels(Kseis(iter), Model(iter));
                    fig_knl = plot_kernels_rho_mu_lambda_relative(K_rel);
                    figname = ['../output/iter',num2str(iter),'.kernels.relative.rho-mu-lambda.png'];
                    titel = [project_name,': seismic kernels (relative rhomulambda) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    print(fig_knl,'-dpng','-r400',figname);
                    close(fig_knl);
                case 'rhovsvp'
                    [Kabs, K_reltemp] = calculate_other_kernels(Kseis(iter), Model(iter));
                    %                 K_reltemp = change_parametrisation_kernels('rhomulambda','rhovsvp',K_rel, Model(i));
                    fig_knl = plot_kernels_rho_vs_vp_relative(K_reltemp);
                    titel = [project_name,': seismic kernels (relative rhovsvp) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.relative.rho-vs-vp.png'];
                    print(fig_knl,'-dpng','-r400',figname);
                    close(fig_knl);
                    fig_knl = plot_kernels_rho_vs_vp(Kabs);
                    titel = [project_name,': seismic kernels (absolute rhovsvp) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.abs.rho-vs-vp.png'];
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
    if (iter<InvProps.niter)
%    if i>1

        if strcmp(use_grav,'yes')
            % determine weight of relative kernels
            w_Kseis = 1;
            w_Kg = 1e7;
            
            disp ' ';
            disp '---'
            disp(['seismic kernel 98th prctile:        ',num2str(prctile(abs(Kseis(iter).rho.total(:)),98))]);
            disp(['gravity kernel 98th prctile:        ',num2str(prctile(abs(Kg{iter}(:)),98))]);
            disp '---'
            disp(['norm seismic kernel:                ',num2str(norm(Kseis(iter).rho.total(:)),'%3.2e')]);
            disp(['norm gravity kernel:                ',num2str(norm(Kg{iter}(:)),'%3.2e')]);
            disp '---'
            InvProps.verhouding98th(iter) = prctile(abs(Kseis(iter).rho.total(:)),98) / prctile(abs(Kg{iter}(:)),98);
            InvProps.verhouding(iter) = norm(Kseis(iter).rho.total(:)) / norm(Kg{iter}(:));
            disp(['the ratio of seis and grav kernel 98th prctiles: ',num2str(InvProps.verhouding98th(iter),'%3.2e')]);
            disp(['the ratio of seis and grav kernel norms:         ',num2str(InvProps.verhouding(iter),'%3.2e')]);
            disp(['the ratio of grav and seis weights:              ',num2str(w_Kg/w_Kseis,'%3.2e')]);
            disp '---'
            
            % combine seismic and gravity kernels
            disp ' ';
            disp(['iter ',num2str(iter),': combining gravity and seismic kernels']);
            Ktest = change_parametrisation_kernels('rhomulambda',parametrisation,Kseis(iter),Model(iter));
            Ktest.rho2.total = w_Kseis * Ktest.rho2.total  +  w_Kg * Kg{iter};
            Ktest1 = change_parametrisation_kernels(parametrisation,'rhomulambda', Ktest,Model(iter));
            K_total(iter).rho.total    = Ktest1.rho.total;
            K_total(iter).mu.total     = Ktest1.mu.total;
            K_total(iter).lambda.total = Ktest1.lambda.total;
            clearvars('Ktest', 'Ktest1');
        else
            K_total(iter) = Kseis(iter);
        end
%     end

            % empty the big variables so that the computer doesn't slow down.
            clearvars('u_fw');
            clearvars('v_fw');

        % calculate the step length and model update
        disp ' ';
        disp(['iter ',num2str(iter),': calculating step length']);
        if iter==1;
            % basis for new step is initial step length
            [InvProps.step(iter), fig_lnsrch] = calculate_step_length(InvProps.stepInit,iter, ...
                InvProps.misfit(iter), InvProps.misfit_seis, InvProps.misfit_g, ...
                Model(iter),K_total(iter),v_obs, g_obs);
        elseif iter>1;
            % basis for new step length is previous one.
            [InvProps.step(iter), fig_lnsrch] = calculate_step_length(InvProps.step(iter-1),iter, ...
                InvProps.misfit(iter), InvProps.misfit_seis, InvProps.misfit_g, ...
                Model(iter),K_total(iter),v_obs, g_obs);
        end        
        % save linesearch figure
%         figname = ['../output/iter',num2str(i),'.step-linesearch.png'];
%         titel = [project_name,': linesearch for iter ',num2str(i)];
%         mtit(fig_lnsrch,titel, 'xoff', 0.001, 'yoff', 0.00001);
%         print(fig_lnsrch,'-dpng','-r400',figname);
%         close(fig_lnsrch);
        clearvars fig_lnsrc;

        % actual model update
        disp ' ';
        disp(['iter ',num2str(iter),': updating model']);
        Model(iter+1) = update_model(K_total(iter),InvProps.step(iter),Model(iter),parametrisation);

        
        
        %% HARD CONSTRAINTS
        % apply hard constraints
        if(strcmp(apply_hc,'yes'))
            % -> no negative velocities
            % -> mass of the Earth and/or its moment of inertia
            switch parametrisation
                case 'rhomulambda'
                    [Model(iter+1).rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model(iter+1).rho,axrot);
                case 'rhovsvp'
                    Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model(iter+1));
                    [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                    Model(iter+1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
                otherwise
                    error('the parametrisation of the inversion was not recognised')
            end
            figname = ['../output/iter',num2str(iter),'.hard-constraints-rhoupdate.png'];
            titel = [project_name,': hard constraints update for iter ',num2str(iter)];
            mtit(fig_hcupdate,titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_hcupdate,'-dpng','-r400',figname);
            close(fig_hcupdate);
            clearvars fig_rhoupdate
        end
        
%         %% calculating model norm
%         InvProps.modeldifnorm(i) =   norm( (Model(i+1).rho(:) - Model(i).rho(:)) ./ Model(i).rho(:) ) ...
%                           + norm( (Model(i+1).mu(:)  - Model(i).mu(:))  ./ Model(i).mu(:) ) ...
%                           + norm( (Model(i+1).lambda(:) - Model(i).lambda(:)) ./ Model(i).lambda(:) );
        
    end
    
    %% OUTPUT:
    
%     % look at angles between gravity kernels in consecutive iterations
%     InvProps.angleKg(i).kernelskeerelkaar = Kg{i-1}(:)'*Kg{i}(:);
%     InvProps.angleKg(i).norm = norm(Kg{i-1}(:))*norm(Kg{i}(:));
%     InvProps.angleKg(i).angle = acos(angle(i).kernelskeerelkaar / angle(i).norm);
    
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
% save(savename, 'InvProps.misfit', '-v7.3');

% plot nice vector plot of misfit evolution + real, start, end model
plot_models_vector;

% useful output
for iter = 1:InvProps.niter
InvProps.misfitseis(iter) = InvProps.misfit_seis{iter}.normd;
InvProps.misfitgrav(iter) = InvProps.misfit_g(iter).normd;
end
% for iter = 1:InvProps.niter
% InvProps.misfitgrav(iter) = InvProps.misfit_g(iter).normd;
% end

% % plot+save real, start and end models as vector plots (also inversion result)
% plot_models_vector;

% inversion results with inversion landscape plot
fig_inv2 = plot_inversion_development_landscapeshape(InvProps);
figname = ['../output/inversion_development.',project_name,'.misfit-landscape.png'];
print(fig_inv2,'-dpng','-r400',figname);
figname = ['../output/inversion_development.',project_name,'.misfit-landscape.eps'];
print(fig_inv2,'-deps','-r400',figname);



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