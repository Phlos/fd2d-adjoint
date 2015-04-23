
% preparation
path(path,'../input');
path(path,'../code');
path(path,'../code/propagation');
path(path,'../tools');
path(path,'../tools/misfits');
path(path,'../quivers');
path(path,'../mtit');

% number of iterations
InvProps.niter = 60;
istart = 1;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, axrot, apply_hc, use_grav, fix_velocities, ...
    use_matfile_startingmodel, starting_model, ...
    true_model_type, f_maxlist, change_src_every, ...
    parametrisation, param_plot, rec_g, X, Z, misfit_type, ...
    normalise_misfits, InvProps.stepInit] = get_input_info;


%% welcome

disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp(['INVERSION RUN ', project_name]);
disp(['-- parametrisation:  ', parametrisation]);
disp '-- using seis? ..... YES!!'
if strcmp(use_grav , 'yes')
    disp '-- using grav? ..... YES!!'
else
    disp '-- using grav? ..... no'
end
if strcmp(apply_hc , 'yes')
    disp '-- using h.c.? ..... YES!!'
else
    disp '-- using h.c.? ..... no'
end
if strcmp(fix_velocities , 'yes')
    disp '-- fixing vels? .... YES!!'
else
    disp '-- fixing vels? .... no'
end
disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';


%% OBS

% freq consists of sources.v_obs and sources.frequency
if ~exist('sources','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var')
    if istart == 1
        disp 'no OBS present, preparing obs...';
        [Model_real, sources, t_obs, props_obs, g_obs] = prepare_obs(true_model_type);
    else
        error('There are no observed properties!')
    end
else
    disp 'obs properties all present... proceeding...'
end

if istart == 1
    disp 'saving obs variables to matfile...'
    savename = ['../output/',project_name,'.obs.all-vars.mat'];
    save(savename, 'sources', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');
end

%=========================================================================

%% Inversion preparation

% save initial model (from input_parameters)
Model_start = update_model();

% set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));



% plot initial model
if istart == 1
    fig_mod = plot_model(Model_start,middle,param_plot);
    set(fig_mod,'Renderer','painters');
    titel = [project_name,': model of iter 0'];
    mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
    figname = ['../output/iter0.model.',param_plot,'.png'];
    print(fig_mod,'-dpng','-r400',figname);
    close(fig_mod);
    
    
    % if 1st iter model @ matfile, load matfile
    if strcmp(use_matfile_startingmodel,'yes')
        load(starting_model)
        Model(1) = Model_out;
        clearvars Model_out;
    else
        Model(1) = Model_start;
    end
    
    % apply hard constraints to initial model (if applicable)
    if(strcmp(apply_hc,'yes'))
        % -> no negative velocities
        % -> mass of the Earth and/or its moment of inertia
        %     [Model(1).rho, fig_hcupdate,~,~] = ...
        %                 apply_hard_constraints(props_obs, Model_start.rho,axrot);
        param_applyhc = 'rhovsvp';
        switch param_applyhc;
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
        
    end
    
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
        
        %% DIFF SOURCES USED 

%         change_src_every = 10;       % how many iterations with the same src?
        cse = change_src_every;
        which_src = floor((iter-1)/cse)+1;
        if which_src > length(sources)
            which_src = length(sources);
        end
        freqmax(iter) = sources(which_src).frequency;
        vobs = sources(which_src).v_obs;
        stf{iter} = sources(which_src).stf;

%         if iter <= length(sources);
%             freqmax(iter) = sources(iter).frequency;
%             vobs = sources(iter).v_obs;
%             stf{iter} = sources(iter).stf;
%         else
%             freqmax(iter) = f_maxlist(end);
%             vobs = sources(end).v_obs;
%             stf{iter} = sources(end).stf;
%         end
        

        
        %% MODEL
        
        % plot model
%         fig_mod = plot_model(Model(iter),middle,param_plot);
%         titel = [project_name,': model of iter ', num2str(iter)];
%         mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
%         figname = ['../output/iter',num2str(iter),'.model.',param_plot,'.png'];
%         print(fig_mod,'-dpng','-r400',figname);
%         close(fig_mod);
%         clearvars('fig_mod');
        fig_mod = plot_model_diff(Model(iter),Model_start,param_plot);
        titel = [project_name,': model diff of iter ', num2str(iter), ' and starting model'];
        mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = ['../output/iter',num2str(iter),'.model-diff.',param_plot,'.png'];
        print(fig_mod,'-dpng','-r400',figname);
        close(fig_mod);
        clearvars('fig_mod');
        
        %% misfit
        
        % gravity misfit
%         if strcmp(use_grav,'yes')
            % gravity field of current model
            [g(iter), fig_grav] = calculate_gravity_field(Model(iter).rho, rec_g);
%             figname = ['../output/iter',num2str(i),'.gravity_recordings.png'];
%             titel = [project_name,': gravity field of ', num2str(i), 'th model'];
%             mtit(fig_grav, titel, 'xoff', 0.001, 'yoff', 0.00001);
%             print(fig_grav, '-dpng', '-r400', figname);
            close(fig_grav);
            % comparison to real model:
            if (strcmp(use_grav,'yes'))
            fig_grav_comp = plot_gravity_quivers(rec_g, g(iter), g_obs, X, Z, Model(iter).rho);
            figname = ['../output/iter',num2str(iter),'.gravity_difference.png'];
            titel = [project_name,': gravity diff iter ', num2str(iter), ' - real model'];
            mtit(fig_grav_comp, titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_grav_comp, '-dpng', '-r400', figname);
            close(fig_grav_comp);
            clearvars('fig_mod');
            end
            
            %- calculate gravity misfit:
            [g_src, InvProps.misfit_g(iter)] = make_gravity_sources(g(iter), g_obs);
            
            % NEW AS OF 17-3-2015 - changed as of 25-3-2015
            InvProps.misfit_g(iter).normd = ...
                norm_misfit(InvProps.misfit_g(iter).total, ...
                            normalise_misfits, InvProps.misfit_g(1).total, ...
                            g_obs);
            clearvars norm_misf;
            
            % output for gravity misfit
            sumgobs = sum(g_obs.mag(:) .^2);
            div_by_gobs = InvProps.misfit_g(iter).total / sumgobs;
            disp ' ';
            disp(['gravity misfit for iter ',num2str(iter),': ', ...
                num2str(InvProps.misfit_g(iter).total,'%3.2e')])
            disp(['   fraction of g_obs: ', ...
                num2str(div_by_gobs,'%3.2e')])
            disp(['   normalised ',normalise_misfits, ': ', ...
                num2str(InvProps.misfit_g(iter).normd,'%3.2e')])
            disp ' ';


        

        
        % seismic misfit
        
        % run forward wave propagation
        disp ' ';
        disp(['iter ',num2str(iter),', src ',num2str(which_src),': calculating forward wave propagation']);
%         clearvars u_fw v_fw;

        % NEW as of 24-3-2015 (that stf{iter} is now supplied)
        [v_iter{iter},t,u_fw,v_fw,rec_x,rec_z]=run_forward(Model(iter), stf{iter});
        close(clf);
        close(clf);
        close(clf);
        
        % plot seismogram difference

        fig_seisdif = plot_seismogram_difference(vobs,v_iter{iter},t);
        titel = [project_name,'difference between seismograms iter ', num2str(iter), ' and obs (source ',num2str(which_src),')'];
        mtit(fig_seisdif, titel, 'xoff', 0.001, 'yoff', 0.02);
        figname = ['../output/iter',num2str(iter),'.seisdif.png'];
        print(fig_seisdif,'-dpng','-r400',figname);
        close(fig_seisdif);
        
        % make adjoint sources
        cd ../tools
        disp ' '; disp(['iter ',num2str(iter),': calculating adjoint stf']);
        
        [adstf{iter}, InvProps.misfit_seis{iter}] = calc_misfitseis_adstf(misfit_type,t,v_iter{iter},vobs);
        % NEW AS OF 17-3-2015
            InvProps.misfit_seis{iter}.normd = ...
                norm_misfit(InvProps.misfit_seis{iter}.total, ...
                            normalise_misfits, ...
                            InvProps.misfit_seis{1}.total, vobs);
        
        
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
        for ii = 1:length(vobs)
            comp = fieldnames(vobs{ii});
            for icomp = 1:length(comp)
            sumvobs = sumvobs + sum(vobs{ii}.(comp{icomp}) .^2);
            end
        end
        div_by_vobs = InvProps.misfit_seis{iter}.total / sumvobs;
        disp ' ';
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
                Kg{iter} = norm_kernel(Kg_temp, normalise_misfits, ...
                    InvProps.misfit_g(1).total, g_obs);
                clearvars Kg_temp;
                
                
                %  plot gravity kernel
                figname = ['../output/iter',num2str(iter),'.kernel_grav.rho.png'];
                titel = [project_name, '- Gravity kernel for iter ',num2str(iter)];
                mtit(fig_Kg,titel, 'xoff', 0.001, 'yoff', 0.00001);
                print(fig_Kg,'-dpng','-r400',figname);
                close(fig_Kg);
                clearvars fig_Kg;
            end
        end
        
        
        %% SEISMIC KERNEL
        
        if(iter < InvProps.niter) % kernels only to be calculated when a next iteration will take place.

            
            % run adjoint to obtain seismic kernels
            disp ' ';
            disp(['iter ',num2str(iter),': calculating adjoint wave propagation']);
            cd ../code/
            Kseis_temp = run_adjoint(u_fw,v_fw,adstf{iter},Model(iter));
            
            % normalise kernels
            Kseis(iter) = norm_kernel(Kseis_temp, normalise_misfits, ...
                InvProps.misfit_seis{1}.total, vobs);
%             clearvars Kseis_temp;
            
            
            % plot the kernels
            disp ' ';
            disp(['iter ',num2str(iter),': plotting kernels']);
            cd ../tools/
            
%             [Kabs, K_reltemp] = calculate_other_kernels(Kseis(iter), Model(iter));
            switch parametrisation
                case 'rhomulambda'

%                     % relative rho-mu-lambda
% %                     fig_knl = plot_kernels_rho_mu_lambda_relative(K_reltemp);
%                     fig_knl = plot_kernels(K_reltemp, 'rhomulambda',Model(iter), 'total', 'same', 99.8); 
%                     figname = ['../output/iter',num2str(iter),'.kernels.relative.rho-mu-lambda.png'];
%                     titel = [project_name,' - seismic kernels (relative rhomulambda) for iter ',num2str(iter)];
%                     mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
%                     print(fig_knl,'-dpng','-r400',figname); close(fig_knl);

                    % absolute rho-mu-lambda
%                     fig_knl = plot_kernels_rho_vs_vp_relative(K_reltemp);
                    fig_knl = plot_kernels(Kseis(iter), 'rhomulambda',Model(iter), 'total', 'own', 99.95);
                    titel = [project_name,' - seismic kernels (absolute rho-mu-lambda) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.absolute.rho-mu-lambda.png'];
                    print(fig_knl,'-dpng','-r400',figname); close(fig_knl);

%                     % relative rho-vs-vp
% %                     fig_knl = plot_kernels_rho_mu_lambda(Kseis(iter));
%                     fig_knl = plot_kernels(K_reltemp, 'rhovsvp',Model(iter), 'total', 'same', 99.8);
%                     figname = ['../output/iter',num2str(iter),'.kernels.relative.rho-vs-vp.png'];
%                     titel = [project_name,' - seismic kernels (relative rho-vs-vp) for iter ',num2str(iter)];
%                     mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
%                     print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
                    
                    % absolute rho-vs-vp
%                     fig_knl = plot_kernels_rho_vs_vp(Kabs);
                    fig_knl = plot_kernels(Kseis(iter), 'rhovsvp',Model(iter), 'total', 'own', 99.95);
                    titel = [project_name,' - seismic kernels (absolute rho-vs-vp) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.abs.rho-vs-vp.png'];
                    print(fig_knl,'-dpng','-r400',figname); close(fig_knl);

                case 'rhovsvp'
%                     % relative rho-vs-vp
% %                     fig_knl = plot_kernels_rho_vs_vp_relative(K_reltemp);
%                     fig_knl = plot_kernels(K_reltemp, 'rhovsvp',Model(iter), 'total', 'same', 99.9);
%                     titel = [project_name, ' - seismic kernels (relative rhovsvp) for iter ',num2str(iter)];
%                     mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
%                     figname = ['../output/iter',num2str(iter),'.kernels.relative.rho-vs-vp.png'];
%                     print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
                    % absolute rho-mu-lambda
%                     fig_knl = plot_kernels_rho_vs_vp_relative(K_reltemp);
                    fig_knl = plot_kernels(Kseis(iter), 'rhomulambda',Model(iter), 'total', 'own', 99.95);
                    titel = [project_name,' - seismic kernels (absolute rho-mu-lambda) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.absolute.rho-mu-lambda.png'];
                    print(fig_knl,'-dpng','-r400',figname); close(fig_knl);

                    % absolute rho-vs-vp
%                     fig_knl = plot_kernels_rho_vs_vp(Kabs);
                    fig_knl = plot_kernels(Kseis(iter), 'rhovsvp',Model(iter), 'total', 'own', 99.95);
                    titel = [project_name,' - seismic kernels (absolute rho-vs-vp) for iter ',num2str(iter)];
                    mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
                    figname = ['../output/iter',num2str(iter),'.kernels.abs.rho-vs-vp.png'];
                    print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
                otherwise
                    error('unrecognised parametrisation for kernel plot');
            end
        end
%         clearvars K_reltemp fig_knl;
        
    
    
    %% COMBINE KERNELS
    
    
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
            disp ' '; disp(['iter ',num2str(iter),': combining gravity and seismic kernels']);
            param_addknls = 'rhovsvp';  % gravity kernel doesn't say anything about vs and vp!!!
                                        % should ALWAYS be applied in
                                        % rhovsvp
            Ktest = change_parametrisation_kernels('rhomulambda',param_addknls,Kseis(iter),Model(iter));
            switch param_addknls
                case 'rhomulambda'
                    Ktest.rho.total = w_Kseis * Ktest.rho.total + w_Kg * Kg{iter};
                case 'rhovsvp'
                    Ktest.rho2.total = w_Kseis * Ktest.rho2.total  +  w_Kg * Kg{iter};
                otherwise
                    error('the parametrisation in which kernels are added was unknown');
            end
            Ktest1 = change_parametrisation_kernels(param_addknls,'rhomulambda', Ktest,Model(iter));
            K_total(iter).rho.total    = Ktest1.rho.total;
            K_total(iter).mu.total     = Ktest1.mu.total;
            K_total(iter).lambda.total = Ktest1.lambda.total;
            clearvars('Ktest', 'Ktest1');
        else
            K_total(iter) = Kseis(iter);
        end

    % empty the big variables so that the computer doesn't slow down.
    clearvars u_fw v_fw;

    %% CALC STEP LN and UPDATE MODEL
        % calculate the step length and model update
        disp ' ';  disp(['iter ',num2str(iter),': calculating step length']);
        if iter==1; stapje = InvProps.stepInit;
        elseif iter>1; stapje = InvProps.step(iter-1);
        end
        [InvProps.step(iter), fig_lnsrch,  InvProps.steplnArray{iter}, InvProps.misfitArray{iter} ] = ...
                calculate_step_length(stapje,iter, ...
                InvProps.misfit(iter), InvProps.misfit_seis, InvProps.misfit_g, ...
                Model(iter),K_total(iter),vobs, g_obs, stf{iter}, Model_start);
        clearvars stapje;     

        % save linesearch figure
        figname = ['../output/iter',num2str(iter),'.step-linesearch.png'];
        titel = [project_name,': linesearch for iter ',num2str(iter)];
        mtit(fig_lnsrch,titel)%, 'xoff', 0.001, 'yoff', 0.00001);
        print(fig_lnsrch,'-dpng','-r400',figname);
        close(fig_lnsrch);
        clearvars fig_lnsrc;

        % initial model update (pre hc, pre fix velocities)
        disp ' ';
        disp(['iter ',num2str(iter),': updating model']);
        Model_prehc(iter+1) = update_model(K_total(iter),InvProps.step(iter),Model(iter),parametrisation);

        
        
        %% HARD CONSTRAINTS
        % apply hard constraints
        if(strcmp(apply_hc,'yes'))
            % -> no negative velocities
            % -> mass of the Earth and/or its moment of inertia
            param_applyhc = 'rhovsvp';
            switch param_applyhc
                case 'rhomulambda'
                    [Model(iter+1).rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_prehc(iter+1).rho,axrot);
                case 'rhovsvp'
                    Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model_prehc(iter+1));
                    [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                    Model_prevfix = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
                otherwise
                    error('the parametrisation of the inversion was not recognised')
            end
            figname = ['../output/iter',num2str(iter),'.hard-constraints-rhoupdate.png'];
            titel = [project_name,': hard constraints update for iter ',num2str(iter)];
            mtit(fig_hcupdate,titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_hcupdate,'-dpng','-r400',figname);
            close(fig_hcupdate);
            clearvars fig_rhoupdate
        else
            Model_prevfix = Model_prehc(iter+1);
        end
        
        if(strcmp(fix_velocities,'yes'))
            Model(iter+1) = fix_vs_vp(Model_prevfix, Model_start);
        else
            Model(iter+1) = Model_prevfix;
        end
%         clearvars Model_prevfix Model_prehc;

    end
    
    %% OUTPUT:
    

    % useful output
    if strcmp(use_grav,'no')
        Kg{iter}=NaN;
    end
    InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model);

    if (iter > 1)
        % inversion results with inversion landscape plot
        fig_inv2 = plot_inversion_development_landscapeshape(InvProps, iter);
        figname = ['../output/inversion_development.',project_name,'.misfit-landscape.png'];
        print(fig_inv2,'-dpng','-r400',figname);
        figname = ['../output/inversion_development.',project_name,'.misfit-landscape.eps'];
        print(fig_inv2,'-depsc','-r400',figname);
        close(fig_inv2)
        
        fig_invres = plot_inversion_result(InvProps, iter);
        titel = [project_name,': inversion results'];
        mtit(fig_invres,titel, 'xoff', 0.0000001, 'yoff', 0.03);
        figname = ['../output/inversion_result.',project_name,'.png'];
        print(fig_invres,'-dpng','-r400',figname);
        figname = ['../output/inversion_result.',project_name,'.eps'];
        print(fig_invres,'-depsc','-r400',figname);
        close(fig_invres)
    end
    
    
    
    %% safety
    if mod(iter,10) == 0
        % saving current variables to file (crash safeguard)
        disp 'saving all current variables...'
        close all;
        clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
            'filenm_old', 'filenm_new', 'fig_knl');
        %     exclude_vars = {'u_fw', 'v_fw'};
        savename = ['../output/',project_name,'.all-vars.mat'];
        save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');
    end
    
end







%% WRAP-UP: misfit evolution & saving all variables to file

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|         ...FINISHING UP...         |';
% disp '======================================';

fig_end = plot_model_diff(Model(niter), Model_start, 'rhovsvp');
clim_rounded = 10*round(fig_end.Children(6).CLim(2) / 10);
for ii = 2:2:6; 
%     clim_rounded = 10*round(fig_end.Children(6).CLim(2) / 10);
    fig_end.Children(ii).CLim = [-clim_rounded clim_rounded]; 
end
titel = [project_name, ' - Final - starting model'];
mtit(fig_end, titel); %, 'xoff', 0.001, 'yoff', 0.02);
figname = ['../output/out.model-diff-final-starting.',param_plot,'.png'];
print(fig_end,'-dpng','-r600',figname);
close(fig_end);


% % plot nice vector plot of misfit evolution + real, start, end model
% disp 'plotting nice vector figures of models and misfit evo'
% plot_models_vector;



% disp 'saving all current variables...'
% clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
%     'filenm_old', 'filenm_new', 'fig_knl');
% savename = ['../output/',project_name,'.all-vars.mat'];
% save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');

% disp ' ';
% disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '|               DONE!                |'
disp '======================================';
disp ' ';
disp(['  (this was inversion ',project_name,')']);