
% preparation
path(path,'../input');
path(path,'../code');
path(path,'../code/propagation');
path(path,'../tools');
path(path,'../tools/misfits');
path(path,'../quivers');
path(path,'../mtit');

% number of iterations
InvProps.niter = 10;
istart = 3;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, axrot, apply_hc, use_grav, parametrisation, rec_g, X, Z, misfit_type, normalise_misfits, InvProps.stepInit] = get_input_info;



%% OBS

% MAKE SURE V_OBS IS PRESENT!!! AND SAVED!!!
%[Model_real, v_obs, t_obs, props_obs, g_obs] = prepare_obs(modelnr);
 


% saving the observed variables
if istart == 1
    disp 'saving obs variables to matfile...'
    savename = ['../output/',project_name,'.obs.all-vars.mat'];
    save(savename, 'v_obs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');
else
    if ~exist('v_obs','var') || ~exist('t_obs','var') || ...
            ~exist('Model_real','var') || ~exist('props_obs','var') || ...
            ~exist('g_obs','var')
        error('There are no observed properties!')
    end
end

%% Inversion preparation

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



% save initial model (from input_parameters)
Model_start = update_model();

% set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));

% plot initial model
fig_mod = plot_model(Model_start,middle,parametrisation);
set(fig_mod,'Renderer','painters');
titel = [project_name,': model of iter 0'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter0.model.',parametrisation,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
close(fig_mod);

% apply hard constraints to initial model (if applicable)
if(strcmp(apply_hc,'yes'))
    % -> no negative velocities
    % -> mass of the Earth and/or its moment of inertia
    Model(1) = Model_start;
    [Model(1).rho, fig_hcupdate,~,~] = ...
                apply_hard_constraints(props_obs, Model_start.rho,axrot);
%     switch parametrisation
%         case 'rhomulambda'
%             [Model(1).rho, fig_hcupdate,~,~] = ...
%                 apply_hard_constraints(props_obs, Model_start.rho,axrot);
%         case 'rhovsvp'
%             Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model_start);
%             [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
%                 apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
%             Model(1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
%         otherwise
%             error('the parametrisation of the inversion was not recognised')
%     end
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

%% start of iterations

for iter = istart : InvProps.niter;
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
        
        fig_mod = plot_model_diff(Model(iter),Model_start,parametrisation);
        titel = [project_name,': model diff of iter ', num2str(iter), ' and starting model'];
        mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = ['../output/iter',num2str(iter),'.model-diff.',parametrisation,'.png'];
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
            [g_src, InvProps.misfit_g(iter)] = make_gravity_sources(g(iter), g_obs);
            
            % NEW AS OF 17-3-2015
            InvProps.misfit_g(iter).normd = ...
                norm_misfit(InvProps.misfit_g(iter).total, ...
                            InvProps.misfit_g(1).total, normalise_misfits);
            clearvars norm_misf;
            
            % output for gravity misfit
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = InvProps.misfit_g(iter).total / sumgobs;
            disp ' ';
            disp(['gravity misfit for iter ',num2str(iter),': ', ...
                num2str(InvProps.misfit_g(iter).total,'%3.2e')])
            disp(['   fraction of g_obs ',num2str(iter),': ', ...
                num2str(div_by_gobs,'%3.2e')])
            disp(['   normalised ',normalise_misfits, ': ', ...
                num2str(InvProps.misfit_g(iter).normd,'%3.2e')])
            disp ' ';


        

        
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
        disp ' '; disp(['iter ',num2str(iter),': calculating adjoint stf']);
        
        [adstf{iter}, InvProps.misfit_seis{iter}] = calc_misfitseis_adstf(misfit_type,t,v_iter{iter},v_obs);
        % NEW AS OF 17-3-2015
            InvProps.misfit_seis{iter}.normd = ...
                norm_misfit(InvProps.misfit_seis{iter}.total, ...
                            InvProps.misfit_seis{1}.total, normalise_misfits);
        
        
%         % plot stf to adstf for one station
%         stationnr = 5;
%         fig_stftoadstf = plot_seismogram_difference(stf, disp, vel, adstf, t);
%         titel = [project_name,': stf to adstf at station', num2str(stationnr), ', iter ', num2str(i), ' and obs'];
%         mtit(fig_stftoadstf, titel, 'xoff', 0.001, 'yoff', 0.02);
%         figname = ['../output/iter',num2str(i),'.stf-to-adstf-station',num2str(stationnr),'.png'];
%         print(fig_stftoadstf,'-dpng','-r400',figname);
%         close(fig_stftoadstf);

       
        %% output total misfit
        
        disp ' ';
        disp '=========================================='
        disp(['           misfit ITER ',num2str(iter)]);
        
        % gravity misfit
        if strcmp(use_grav,'yes')
            sumgobs = sum(g_obs.x(:) .^2) + sum(g_obs.z(:) .^2);
            div_by_gobs = InvProps.misfit_g(iter).total / sumgobs;
            disp(['GRAVITY misfit FOR ITER ',num2str(iter),': ', ...
                num2str(InvProps.misfit_g(iter).total,'%3.2e')])
            disp(['   fraction of g_obs:       ', ...
                num2str(div_by_gobs,'%3.2e')])
            disp(['   normalised ',normalise_misfits,':  ', ...
                num2str(InvProps.misfit_g(iter).normd,'%3.2e')])
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
%         disp ' ';
        disp(['SEISMIC misfit FOR ITER ',num2str(iter),': ', ...
            num2str(InvProps.misfit_seis{iter}.total,'%3.2e')])
        disp(['   fraction of v_obs:       ', ...
            num2str(div_by_vobs,'%3.2e')])
%         if strcmp(normalise_misfits,'byfirstmisfit')
        disp(['   normalised ',normalise_misfits,':  ', ...
            num2str(InvProps.misfit_seis{iter}.normd,'%3.2e')])
%         end
        disp ' ';
        
        if strcmp(use_grav,'yes')
            InvProps.misfit(iter) = InvProps.misfit_seis{iter}.normd + InvProps.misfit_g(iter).normd;
        else
            InvProps.misfit(iter) = InvProps.misfit_seis{iter}.normd;
        end
        
        disp(['TOTAL misfit FOR ITER ',num2str(iter),': ', ...
            num2str(InvProps.misfit(iter),'%3.2e')])
%         disp ' ';
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

                % calculating the gravity kernel
                if iter == 1
                    [Kg_temp, fig_Kg] = compute_kernels_gravity(g_src,rec_g,'no'); % 'no' is for plotting gravity kernel update
                else
                    [Kg_temp, fig_Kg] = compute_kernels_gravity(g_src,rec_g,'no'); % 'no' is for plotting gravity kernel update
                end
                
                % normalising the gravity kernel
                Kg{iter} = norm_kernel(Kg_temp, InvProps.misfit_g(1).total, normalise_misfits);
                clearvars Kg_temp;
                
%                 normKg = InvProps
%                 if strcmp(normalise_misfits,'byfirstmisfit')
%                     normKg = 1.0 / InvProps.misfit_g(1).total;
%                 else
%                     normKg = 1.0;
%                 end
                
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
        
        if(iter < InvProps.niter) % kernels only to be calculated when a next iteration will take place.
            %  if (i>1)
            % normalising the seismic kernel
%             if strcmp(normalise_misfits,'byfirstmisfit')
%                 InvProps.normKseis = 1.0 / InvProps.misfit_seis{1}.total;
%             else
%                 normKseis = 1.0;
%             end
            
            % run adjoint to obtain seismic kernels
            disp ' ';
            disp(['iter ',num2str(iter),': calculating adjoint wave propagation']);
            cd ../code/
            Kseis_temp = run_adjoint(u_fw,v_fw,adstf{iter},Model(iter));
            
            % normalise kernels
            Kseis(iter) = norm_kernel(Kseis_temp, InvProps.misfit_seis{1}.total, normalise_misfits);
            clearvars Kseis_temp;
            
            
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
        
    
    
    %% COMBINE KERNELS & UPDATE MODEL
    
    
    % only update the model if we're going to a next model
    if (iter<InvProps.niter)
%    if i>1

        if strcmp(use_grav,'yes')
            % determine weight of relative kernels
            w_Kseis = 1;
            w_Kg = 1;
%             w_Kg = 1e7; % was necessary before discovery KMP problem: spatial delta
            
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
%             disp(['the ratio of grav and seis weights:              ',num2str(w_Kg/w_Kseis,'%3.2e')]);
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


    % empty the big variables so that the computer doesn't slow down.
    clearvars u_fw v_fw;

        % calculate the step length and model update
        disp ' ';  disp(['iter ',num2str(iter),': calculating step length']);
        if iter==1; stapje = InvProps.stepInit;
        elseif iter>1; stapje = InvProps.step(iter-1);
        end
        [InvProps.step(iter), fig_lnsrch,  InvProps.steplnArray{iter}, InvProps.misfitArray{iter} ] = ...
                calculate_step_length(stapje,iter, ...
                InvProps.misfit(iter), InvProps.misfit_seis, InvProps.misfit_g, ...
                Model(iter),K_total(iter),v_obs, g_obs);
        clearvars stapje;     

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
        

    end
    
    %% OUTPUT:
    

    % useful output
    InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model);

    if (iter > 1)
    % inversion results with inversion landscape plot
    fig_inv2 = plot_inversion_development_landscapeshape(InvProps, iter);
    figname = ['../output/inversion_development.',project_name,'.misfit-landscape.png'];
    print(fig_inv2,'-dpng','-r400',figname);
    figname = ['../output/inversion_development.',project_name,'.misfit-landscape.eps'];
    print(fig_inv2,'-deps','-r400',figname);
    close(fig_inv2)
    
    fig_invres = plot_inversion_result(InvProps, iter);
    figname = ['../output/inversion_result.',project_name,'.png'];
    print(fig_invres,'-dpng','-r400',figname);
    figname = ['../output/inversion_result.',project_name,'.eps'];
    print(fig_invres,'-deps','-r400',figname);
    close(fig_invres)
    
    end
    
    
    
    %% safety
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


% % plot nice vector plot of misfit evolution + real, start, end model
% disp 'plotting nice vector figures of models and misfit evo'
% plot_models_vector;



% disp 'saving all current variables...'
% clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
%     'filenm_old', 'filenm_new', 'fig_knl');
% savename = ['../output/',project_name,'.all-vars.mat'];
% save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|               DONE!                |'
disp '======================================';
disp ' ';
disp(['  (this was inversion ',project_name,')']);