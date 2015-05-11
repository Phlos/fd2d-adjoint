function [misfit_total, misfit_seis, misfit_grav] = calc_misfits(Model, ...
    g_obs, misfit_grav_initial,       stf, v_obs, misfit_seis_initial)

% calculates the normalised misfit for a given model: total, seis and grav
% -- Nienke Blom, 29-4-2015

    %% preparation
    input_parameters;

    %% gravity misfit
    
    if strcmp(use_grav,'yes')
        [g_try, fig_grav] = calculate_gravity_field(Model.rho, rec_g);
        close(fig_grav);

        %- calculate gravity misfit:
        [~, misf_g_test] = make_gravity_sources(g_try, g_obs);
    %    disp(['misfit_g_test.total ', misf_g_test.total]);
        misfit_grav = norm_misfit(misf_g_test.total, normalise_misfits, ...
            misfit_grav_initial, g_obs);
    else
        misfit_grav = NaN;
    end
    
    %% seismic misfit
    
    %- run forward update in given model
    [v_try,t,~,~,~,~] = run_forward(Model, stf);
    
%     % plot seismogram difference:
%     fig_seisdif = plot_seismogram_difference(v_obs, v_try, t);
%     figname = ['../output/iter.current.calcstepln.seisdif.',num2str(steptry),'.png'];
%     print(fig_seisdif,'-dpng','-r400',figname);
%     close(fig_seisdif);
    
    %- calculate seismic misfit: 
    [~, misf_s_test] = calc_misfitseis_adstf(misfit_type,t,v_try,v_obs);
    misfit_seis = norm_misfit(misf_s_test.total, normalise_misfits, ...
                                misfit_seis_initial, v_obs);



    %% combine the misfits
    if strcmp(use_grav,'yes')
        misfit_total = misfit_seis + misfit_grav;
    else
        misfit_total =  misfit_seis;
    end
    
    if strcmp(use_grav,'yes')
        disp(['seismic misfit:  ', num2str(misfit_seis,'%3.2e')])
        disp(['gravity misfit:  ', num2str(misfit_grav,'%3.2e')])
        disp(['total misfit:    ', num2str(misfit_total,'%3.2e')])
    else
        disp(['total (=seis) misfit: ', num2str(misfit_total,'%3.2e')]);
    end

end