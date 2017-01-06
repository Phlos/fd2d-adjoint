function [misfit_init, misfit_init_traces]  = calc_initial_misfits(Model_init, sObsPerFreq, g_obs)

input_parameters;
if strcmp(use_seis, 'yesseis')
for ifrq = 1:length(sObsPerFreq)
    
    sEventInfo = sObsPerFreq(ifrq).sEventInfo; % contains src locations etc
    sEventObs = sObsPerFreq(ifrq).sEventObs;  % contains v_obs
    
    disp(['calculating initial model misfit for freq. ',num2str(ifrq)]);
    [misfit_int.total, misfit_int.seis, misfit_int.grav] = ...
        calc_misfits(Model_init, g_obs, 0, sEventInfo, sEventObs, 0, ...
        'noplot','notext');
    misfit_init(ifrq) = misfit_int;
    misfit_init_traces{ifrq} = sEventObs;
    
end

else
    disp 'calculating initial model misfit for gravity only';
    [misfit_init.total, misfit_init.seis, misfit_init.grav] = ...
        calc_misfits(Model_init, g_obs, 0, NaN, NaN, 0, ...
        'noplot','notext');
end

end
