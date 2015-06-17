function misfit_init = calc_initial_misfits(Model_init, sObsPerFreq, g_obs)

% calculates the initial misfits for all supplied frequencies

for ifrq = 1:length(sObsPerFreq)
    
    sEventInfo = sObsPerFreq(ifrq).sEventInfo; % contains src locations etc
    sEventObs = sObsPerFreq(ifrq).sEventObs;  % contains v_obs
    
    disp(['calculating initial model misfit for freq. ',num2str(ifrq)]);
    [misfit_int.total, misfit_int.seis, misfit_int.grav] = ...
        calc_misfits(Model_init, g_obs, 0, sEventInfo, sEventObs, 0, ...
        'noplot','notext');
    misfit_init(ifrq) = misfit_int;
    
end

end