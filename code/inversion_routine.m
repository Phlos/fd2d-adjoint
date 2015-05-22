
%% preparation

% number of iterations
InvProps.niter = 60;
istart = 16;

niter = InvProps.niter;

% obtain useful parameters from input_parameters
[project_name, axrot, apply_hc, use_grav, fix_velocities, ...
    use_matfile_startingmodel, starting_model, bg_model_type,...
    true_model_type, f_maxlist, change_freq_every, ...
    parametrisation, param_plot, rec_g, X, Z, misfit_type, ...
    normalise_misfits, InvProps.stepInit] = get_input_info;


% project folder
output_path = ['../output/',project_name,'/'];
mkdir('../output/',project_name)

% save input_parameters to this path
savename = [output_path,project_name,'.input_parameters.m'];
copyfile('../input/input_parameters.m',savename)

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

% freq consists of freqs.v_obs and freqs.frequency
if ~exist('sObsPerFreq','var') || ~exist('t_obs','var') || ...
        ~exist('Model_real','var') || ~exist('props_obs','var') || ...
        ~exist('g_obs','var')
    if istart == 1
        disp 'no OBS present, preparing obs...';
        [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path,true_model_type);
    else
        error('iter > 1 but there are no observed properties!')
    end
else
    disp 'obs properties all present... proceeding...'
end


savename = [output_path,'/obs.all-vars.mat'];
if ~(exist(savename, 'file'))
    disp 'saving obs variables to matfile...'
    save(savename, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');
end

%=========================================================================

%% Model preparation

disp 'preparing inversion starting model...'

% save initial model (from input_parameters)
Model_start = update_model();

% save background model (for plotting purposes)
Model_bg = update_model(bg_model_type);

% set the background values for plot_model to mode of the initial model
middle.rho    = mode(Model_start.rho(:));
middle.mu     = mode(Model_start.mu(:));
middle.lambda = mode(Model_start.lambda(:));


% plot initial model
if istart == 1
    disp 'plotting initial model'
    fig_mod = plot_model_diff(Model_start,Model_bg,param_plot);
    set(fig_mod,'Renderer','painters');
    titel = [project_name,': starting model'];
    mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
    figname = [output_path,'/iter0.starting-model.',param_plot,'.png'];
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
        param_applyhc = 'rhovsvp';
        switch param_applyhc;
            case 'rhomulambda'
                [Model(1).rho, fig_hcupdate,~,~] = ...
                    apply_hard_constraints(props_obs, Model(1).rho,axrot);
            case 'rhovsvp'
                Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model(1));
                [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                    apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                Model(1) = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
            otherwise
                error('the parametrisation of the inversion was not recognised')
        end
        figname = [output_path,'/iter0.hard-constraints-rhoupdate.png'];
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
        
 
        
        %% MODEL
        
        % plot model
%         fig_mod = plot_model(Model(iter),middle,param_plot);
%         titel = [project_name,': model of iter ', num2str(iter)];
%         mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
%         figname = [output_path,'/iter',num2str(iter),'.model.',param_plot,'.png'];
%         print(fig_mod,'-dpng','-r400',figname);
%         close(fig_mod);
%         clearvars('fig_mod');
        fig_mod = plot_model_diff(Model(iter),Model_bg,param_plot);
        titel = [project_name,': model diff of iter ', num2str(iter), ' and bg model'];
        mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = [output_path,'/iter',num2str(iter,'%03d'),'.model-diff.',param_plot,'.png'];
        print(fig_mod,'-dpng','-r400',figname);
        close(fig_mod);
        clearvars('fig_mod');
        
        %% get current frequency and its source & obs info

        cfe = change_freq_every;
        whichFrq = floor((iter-1)/cfe)+1;
        if whichFrq > length(sObsPerFreq)
            whichFrq = length(sObsPerFreq);
        end
        freqmax(iter) = sObsPerFreq(whichFrq).f_max;
        freqmin(iter) = sObsPerFreq(whichFrq).f_min;
        
        nsrc = length(sObsPerFreq(whichFrq).sEventObs);
        
        % event info = source locations & (unfiltered) stf
        % event obs = observed seismograms per source per receiver
        sEventInfo = sObsPerFreq(whichFrq).sEventInfo; % contains src locations etc
        sEventObs = sObsPerFreq(whichFrq).sEventObs;  % contains v_obs
 
        
        %% misfit
        
        % initial misfit  ---  actually this should only be done if iter>1
        if iter >1
            if(~exist('misfit_init', 'var') || length(misfit_init) < whichFrq || isempty(misfit_init(whichFrq).total) )
                % do somethin with calculate seis misfit for that source @mod1
                disp(['calculating initial model misfit for freq. ',num2str(whichFrq)]);
                [misfit_int.total, misfit_int.seis, misfit_int.grav] = ...
                    calc_misfits(Model(1), g_obs, 0, sEventInfo, sEventObs, 0, ...
                    'noplot','notext');
                misfit_init(whichFrq) = misfit_int;
            end
        end
        if iter ==1;
            misfit_init(whichFrq).total = 0;
            misfit_init(whichFrq).seis = 0;
            misfit_init(whichFrq).grav = 0;
        end
        
        % current misfit
        disp(['calculating current misfit for iter ',num2str(iter)]);

        [misfit_total, misfit_seis, misfit_grav, ...
        giter, g_src, sEventRecIter, sEventAdstfIter] = calc_misfits(Model(iter), ...
                                            g_obs, misfit_init(whichFrq).grav , ...
                                            sEventInfo, sEventObs, misfit_init(whichFrq).seis, ...
                                            'yessavefields','yessaveplots');
        if iter ==1;
            misfit_init(whichFrq).grav = misfit_grav;
            misfit_grav = 1;
            misfit_init(whichFrq).seis = misfit_seis;
            misfit_seis = 1;
            misfit_init(whichFrq).total = misfit_total;
            misfit_total = 2;
        end
        InvProps.misfit(iter) = misfit_total;
        InvProps.misfitseis(iter) = misfit_seis;
        InvProps.misfitgrav(iter) = misfit_grav;
        g(iter) = giter;
        sEventRec{iter} = sEventRecIter;
        sEventAdstf{iter} = sEventAdstfIter;
        
        
        % gravity misfit
        disp ' ';
        disp(['GRAVITY misfit FOR ITER ',num2str(iter,'%2u'),':   ', ...
            num2str(InvProps.misfitgrav(iter),'%3.2e')])
        disp(['   normalised ',normalise_misfits, ...
            ' -- 1st misfit = ',num2str(misfit_init(whichFrq).grav,'%3.2e'),')'])
        disp ' ';
        % seismic misfit
        disp(['SEISMIC misfit FOR ITER ',num2str(iter,'%2u'),':   ', ...
            num2str(InvProps.misfitseis(iter),'%3.2e')])
        disp(['   normalised ',normalise_misfits, ...
            ' -- 1st misfit = ',num2str(misfit_init(whichFrq).seis,'%3.2e'),')'])
        disp ' ';
        % total misfit
        disp(['TOTAL misfit FOR ITER ',num2str(iter,'%2u'),':     ', ...
            num2str(InvProps.misfit(iter),'%3.2e')])
%         disp ' ';
        disp '=========================================='
        


        % plot misfit evolution
        fig_misfit = plot_misfit_evolution(InvProps);
        figname = [output_path,'/misfit-evolution.pdf'];
        mtit(fig_misfit, project_name, 'xoff', 0.001, 'yoff', 0.04);
        print(fig_misfit,'-dpdf','-r400',figname);
        close(fig_misfit);
        clearvars fig_misfit

        
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
                    misfit_init(whichFrq).grav);
                clearvars Kg_temp;
                
                
                %  plot gravity kernel
                figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernel_grav.rho.png'];
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
            disp(['iter ',num2str(iter),': calculating seismic kernels']);
            [Kseis_temp, sEventKnls{iter}] = run_adjoint_persource(Model(iter), sEventAdstf{iter});
            
            % normalise kernels
            Kseis(iter) = norm_kernel(Kseis_temp, normalise_misfits, ...
                misfit_init(whichFrq).seis);
%             clearvars Kseis_temp;
            
            % plot the kernels
            disp ' ';
            disp(['iter ',num2str(iter),': plotting kernels']);
            cd ../tools/
            
            % absolute rho-mu-lambda
            fig_knl = plot_kernels(Kseis(iter), 'rhomulambda',Model(iter), 'total', 'own', 99.95);
            titel = [project_name,' - iter ',num2str(iter), ' seismic kernels (abs rho-mu-lambda)'];
            mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
            figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels.abs.rho-mu-lambda.png'];
            print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
            % absolute rho-vs-vp
            fig_knl = plot_kernels(Kseis(iter), 'rhovsvp',Model(iter), 'total', 'own', 99.95);
            titel = [project_name, ' - iter ',num2str(iter), ' seismic kernels (abs rho-vs-vp)'];
            mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
            figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels.abs.rho-vs-vp.png'];
            print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
            
%             % plot subkernels
%             for isrc = 1:numel(sEventKnls{iter})
%                 Kernel = sEventKnls{iter}(isrc);
%                 % absolute rho-mu-lambda
%                 fig_knl = plot_kernels(Kernel, 'rhomulambda',Model(iter), 'total', 'own', 99.95);
%                 titel = [project_name,' - iter ',num2str(iter),' src ',num2str(isrc),' - seismic SUBkernel (abs rho-mu-lambda)'];
%                 mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
%                 figname = [output_path,'/iter',num2str(iter,'%03d'),'.subkernels-src',num2str(isrc,'%02d'),'.abs.rho-mu-lambda.png'];
%                 print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
%                 % absolute rho-vs-vp
%                 fig_knl = plot_kernels(Kernel, 'rhovsvp',Model(iter), 'total', 'own', 99.95);
%                 titel = [project_name,' - iter ',num2str(iter),' src ',num2str(isrc),' - seismic SUBkernel (absolute rho-vs-vp)'];
%                 mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
%                 figname = [output_path,'/iter',num2str(iter,'%03d'),'.subkernels-src',num2str(isrc,'%02d'),'.abs.rho-vs-vp.png'];
%                 print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
%             end; clearvars Kernel

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
                    Krho.seis = filter_kernels(Ktest.rho.total,5);
                    Krho.grav = filter_kernels(Kg{iter},5);
                    Ktest.rho.total = w_Kseis * Ktest.rho.total + w_Kg * Kg{iter};
                    Krho.together.total = filter_kernels(Ktest.rho.total,5);
                case 'rhovsvp'
                    Krho.seis = filter_kernels(Ktest.rho2.total,5);
                    Krho.grav = filter_kernels(Kg{iter},5);
                    Ktest.rho2.total = w_Kseis * Ktest.rho2.total  +  w_Kg * Kg{iter};
                    Krho.together = filter_kernels(Ktest.rho2.total,5);
                otherwise
                    error('the parametrisation in which kernels are added was unknown');
            end

            
            % saving the total kernel
            Ktest1 = change_parametrisation_kernels(param_addknls,'rhomulambda', Ktest,Model(iter));
            K_total(iter).rho.total    = Ktest1.rho.total;
            K_total(iter).mu.total     = Ktest1.mu.total;
            K_total(iter).lambda.total = Ktest1.lambda.total;
            clearvars('Ktest', 'Ktest1');
        else
            K_total(iter) = Kseis(iter);
        end
        
       % plotting the total kernel in (absolute) rhovsvp
       fig_knl = plot_kernels(K_total(iter), 'rhovsvp',Model(iter), 'total', 'own', 99.95);
       titel = [project_name,' - TOTAL kernels (absolute rho-vs-vp) for iter ',num2str(iter)];
       mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
       figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels-total.abs.rho-vs-vp.png'];
       print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
       
       % plotting the density kernels seis + grav = total
       if strcmp(use_grav, 'yes')
           fig_Krho = plot_model(Krho);
           maks = prctile(abs([Krho.seis(:); Krho.grav(:);Krho.together(:)]),99.5);
           for ii = 2:2:6
               fig_Krho.Children(ii).CLim = [-maks maks];
           end
           titel = [project_name,' - buildup of density kernel - iter ',num2str(iter)];
           mtit(fig_Krho,titel, 'xoff', 0.001, 'yoff', 0.04);
           figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernel-rho-buildup.png'];
           print(fig_Krho,'-dpng','-r400',figname); close(fig_Krho);
       end


    %% CALC STEP LN and UPDATE MODEL
        % calculate the step length and model update
        disp ' ';  disp(['iter ',num2str(iter),': calculating step length']);
        if iter==1; stapje = InvProps.stepInit;
        elseif iter>1; stapje = InvProps.step(iter-1);
        end
        
        % actual step length calculation
        [step, fig_lnsrch, steplnArray, misfitArray] = ...
                calculate_step_length(stapje,iter, ...
                InvProps.misfit(iter), misfit_init(whichFrq), ...
                Model(iter),K_total(iter),g_obs, sEventInfo, sEventObs, Model_start);
        InvProps.step(iter) = step;
        InvProps.steplnArray{iter} = steplnArray;
        InvProps.misfitArray{iter} = misfitArray;
        clearvars stapje step steplnArray misfitArray;

        % save linesearch figure
        figname = [output_path,'/iter',num2str(iter,'%03d'),'.step-linesearch.png'];
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
                    Model_prevfix = Model_prehc(iter+1);
                    [Model_prevfix.rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_prehc(iter+1).rho,axrot);
                case 'rhovsvp'
                    Model_rhovsvp = change_parametrisation('rhomulambda','rhovsvp', Model_prehc(iter+1));
                    [Model_rhovsvp.rho, fig_hcupdate,~,~] = ...
                        apply_hard_constraints(props_obs, Model_rhovsvp.rho,axrot);
                    Model_prevfix = change_parametrisation('rhovsvp','rhomulambda',Model_rhovsvp);
                    clearvars Model_rhovsvp;
                otherwise
                    error('the parametrisation of the inversion was not recognised')
            end
            figname = [output_path,'/iter',num2str(iter,'%03d'),'.hard-constraints-rhoupdate.png'];
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
    if exist('Model_real', 'var')
        InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model, Model_real);
    else
        InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model, Model_real);
    end

    if (iter > 1)
        % inversion results with inversion landscape plot
        fig_inv2 = plot_inversion_development_landscapeshape(InvProps, iter);
        figname = [output_path,'/inversion_development.',project_name,'.misfit-landscape'];
        print(fig_inv2,'-dpng','-r400',[figname,'.png']);
        print(fig_inv2,'-depsc','-r400',[figname,'.eps']);
        close(fig_inv2)
        
        fig_invres = plot_inversion_result(InvProps, iter);
        titel = [project_name,': inversion results'];
        mtit(fig_invres,titel, 'xoff', 0.0000001, 'yoff', 0.03);
        figname = [output_path,'/inversion_result.',project_name];
        print(fig_invres,'-dpng','-r400',[figname,'.png']);
        print(fig_invres,'-depsc','-r400',[figname,'.eps']);
        close(fig_invres)
    end
    
    
    % save all output files to the actual output path
    blips = dir('../output/*iter.current*');
    for ii = 1:numel(blips)
        bestand = blips(ii).name;
        oldfile = ['../output/',bestand];
        newfile = [output_path,strrep(blips(ii).name,'iter.current',['iter',num2str(iter,'%03d')])];
        movefile(oldfile,newfile);
    end
    
    %% safety
    % saving current variables to file (crash safeguard)

%     if mod(iter,10) == 0
        disp 'saving all current variables...'
        close all;
        clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
            'filenm_old', 'filenm_new', 'fig_knl');
        savename = [output_path,'/inversion.all-vars.mat'];
        save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');
%     end
    
end







%% WRAP-UP: misfit evolution & saving all variables to file

disp ' ';
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '======================================';
disp '|         ...FINISHING UP...         |';

fig_end = plot_model_diff(Model(niter), Model_bg, 'rhovsvp');
clim_rounded = 10*round(fig_end.Children(6).CLim(2) / 10);k
for ii = 2:2:6; 
%     clim_rounded = 10*round(fig_end.Children(6).CLim(2) / 10);
    fig_end.Children(ii).CLim = [-clim_rounded clim_rounded]; 
end
titel = [project_name, ' - Final - background model'];
mtit(fig_end, titel); %, 'xoff', 0.001, 'yoff', 0.02);
figname = [output_path,'/out.model-diff-final-bg.',param_plot,'.png'];
print(fig_end,'-dpng','-r600',figname);
close(fig_end);


% % plot nice vector plot of misfit evolution + real, start, end model
% disp 'plotting nice vector figures of models and misfit evo'
% plot_models_vector;


if ~exist([output_path,'/inversion.all-vars.mat'],'file')
    disp '|  (saving all current variables..)  |'
    clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
        'filenm_old', 'filenm_new', 'fig_knl');
    savename = [output_path,'/inversion.all-vars.mat'];
    save(savename, '-regexp', '^(?!(u_fw|v_fw)$).');
end

% disp ' ';
% disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '|               DONE!                |'
disp '======================================';
disp ' ';
disp(['  (this was inversion ',project_name,')']);