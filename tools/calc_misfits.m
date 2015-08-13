function [misfit_total, misfit_seis, misfit_grav, ...
        g_try, g_src, sEventRec, sEventAdstf] = calc_misfits(Model, ...
        g_obs, misfit_grav_initial, sEventInfo, sEventObs, misfit_seis_initial, ...
        varargin)
    
    % calculates the misfit for a given model: total, seis and grav
    % -- Nienke Blom, 29-4-2015
    
    %% preparation
    % plot gravity field or not? Default: yes.
    [plotornot, saveplots, textornot, saveFwdFields] = which_output(varargin(:));
    
    input_parameters;
    [X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
    
    % project folder
    output_path = ['./output/',project_name,'/'];
    
    %     misfit_total = 0;
    nsrc = length(sEventObs);
    
    %% gravity misfit
    
    % calculate gravity field
    g_try = calculate_gravity_field(Model.rho, rec_g, 'noplot');
    
    % plot grav difference
    if  strcmp(plotornot, 'yesplot') % && strcmp(use_grav,'yes')
        fig_grav_comp = plot_gravity_quivers(rec_g, g_try, g_obs, X, Z, Model.rho);
        figname = [output_path,'/iter.current.gravity_difference.png'];
        titel = ['gravity diff of current iter model - real model'];
        if strcmp(saveplots,'yessaveplots')
            mtit(fig_grav_comp, titel, 'xoff', 0.001, 'yoff', 0.00001);
            print(fig_grav_comp, '-dpng', '-r400', figname);
        end
        close(fig_grav_comp);
        clearvars('fig_grav_comp');
    end
    
    %- calculate gravity misfit:
    [g_src, misf_g_test] = make_gravity_sources(g_try, g_obs);
    % disp(['misfit_g_test.total ', misf_g_test.total]);
    misfit_grav = norm_misfit(misf_g_test.total, normalise_misfits, ...
        misfit_grav_initial, g_obs);
    
    % output for gravity misfit
    if strcmp(textornot, 'yestext')
        sumgobs = sum(g_obs.mag(:) .^2);
        div_by_gobs = misf_g_test.total / sumgobs;
        
        disp ' '; disp(['gravity misfit:    ', ...
            num2str(misf_g_test.total,'%3.2e')])
        disp(['   fraction of g_obs:   ', ...
            num2str(div_by_gobs,'%3.2e')])
        disp(['   normalised ',normalise_misfits, ': ', ...
            num2str(misfit_grav,'%3.2e')])
        disp ' ';
    end
    
    %% seismic misfit
    
    if strcmp(use_seis,'yesseis')
        if strcmp(textornot,'yestext')
            disp ' ';
            % disp(['iter ',num2str(iter),', src ',num2str(whichFrq),': calculating forward wave propagation']);
            disp('calculating forward wave propagation');
        end
        
        [sEventRec, fig_seisdif] = run_forward_persource(Model, sEventInfo, sEventObs, saveFwdFields, plotornot);
        
        % save seismograms (+ diff w/ obs) to file
        if strcmp(plotornot,'yesplot')
            for isrc = 1:nsrc
                titel = ['difference between seismograms src ', num2str(isrc), ' and obs'];
                mtit(fig_seisdif(isrc), titel, 'xoff', 0.001, 'yoff', 0.02);
                if strcmp(saveplots,'yessaveplots')
                    figname = [output_path,'/iter.current.seisdif.src',num2str(isrc),'.png'];
                    print(fig_seisdif(isrc),'-dpng','-r400',figname);
                end
                close(fig_seisdif(isrc));
            end
        else
            clearvars fig_seisdif;
        end
        
        %- calculate seismic misfit:
        for isrc = 1:nsrc
            v_rec = sEventRec(isrc).vel;
            v_obs = sEventObs(isrc).vel;
            t = sEventRec(isrc).t;
            
            [sEventAdstf(isrc).adstf, misf_s_test(isrc)] = calc_misfitseis_adstf(misfit_type,t,v_rec,v_obs);
            misfit_s_test(isrc) = norm_misfit(misf_s_test(isrc).total, normalise_misfits, ...
                misfit_seis_initial, v_obs);
            %         % plot stf to adstf for one station
            %         stationnr = 5;
            %         fig_stftoadstf = plot_seismogram_difference(stf, disp, vel, adstf, t);
            %         titel = [project_name,': stf to adstf at station', num2str(stationnr), ', iter ', num2str(i), ' and obs'];
            %         mtit(fig_stftoadstf, titel, 'xoff', 0.001, 'yoff', 0.02);
            %         figname = [output_path,'/iter',num2str(i),'.stf-to-adstf-station',num2str(stationnr),'.png'];
            %         print(fig_stftoadstf,'-dpng','-r400',figname);
            %         close(fig_stftoadstf);
        end
        
        misfit_seis = sum(misfit_s_test);
        
        disp ' '; disp(['seismic misfit (normalised ',normalise_misfits, '):    ', ...
            num2str(misfit_seis,'%3.2e')])
        disp ' ';
    else
        misfit_seis = NaN;
        sEventRec = NaN;
        sEventAdstf = NaN;
        
    end
    
    %% combine the misfits
    misfit_total = 0;
    if strcmp(use_grav,'yes') && strcmp(use_seis, 'yesseis')
        misfit_total = misfit_seis + misfit_grav;
    elseif strcmp(use_seis, 'yesseis')
        misfit_total =  misfit_seis;
    elseif strcmp(use_grav, 'yes')
        misfit_total = misfit_grav;
    end
    
    %     if strcmp(textornot,'yestext')
    if strcmp(use_seis, 'yesseis')
        disp(['seismic misfit:  ', num2str(misfit_seis,'%3.2e')])
    end
    if strcmp(use_grav, 'yes')
        disp(['gravity misfit:  ', num2str(misfit_grav,'%3.2e')])
    end
    disp(['total misfit:    ', num2str(misfit_total,'%3.2e')])
    %     end
    
end

function [plotornot, saveplots, textornot, saveFwdFields] = which_output(args)
    
    % plot figures or not? Output text or not?
    plottext = {'yesplot', 'noplot'};
    saveplottxt = {'yessaveplots', 'nosaveplots'};
    texttext = {'yestext', 'notext'};
    savetext = {'yessavefields', 'nosavefields'};
    
    % default values
    plotornot = 'yesplot';
    saveplots = 'nosaveplots';
    textornot = 'yestext';
    saveFwdFields = 'nosavefields';
    
    for ii = 1:length(args)
        if any(strcmp(args{ii},plottext))
            plotornot = args{ii};
        elseif any(strcmp(args{ii},saveplottxt))
            saveplots = args{ii};
        elseif any(strcmp (args{ii},texttext))
            textornot = args{ii};
        elseif any(strcmp (args{ii},savetext))
            saveFwdFields = args{ii};
        else
            error(['which_output variable input ',num2str(ii),' not recognised']);
        end
    end
    
end