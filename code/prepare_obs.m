function [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path, varargin)

input_parameters;
nfr = length(f_maxlist);


%% Model
Model_real = checkargs(varargin(:));

% Model_real = update_model(modelnr);

% % plotting the real model
% fig_mod_real = plot_model(Model_real, 'rhovsvp');
% mtit(fig_mod_real, 'real model -- rho-vs-vp parametrisation');
% figname = [output_path,'/obs.model.rhovsvp.png'];
% print(fig_mod_real, '-dpng', '-r400', figname);
% close(fig_mod_real);

% real - background model
Mbg = update_model(bg_model_type);
% fig_mod_diff = plot_model_diff(Model_real, Mbg, 'rhovsvp');
% mtit(fig_mod_diff, 'real - bg model -- rho-vs-vp parametrisation');
% figname = [output_path,'/obs.model-real-diff-bg.rhovsvp.png'];
% print(fig_mod_diff, '-dpng', '-r400', figname);
% close(fig_mod_diff);

% h.c. model properties for real model
props_obs = calculate_model_properties(Model_real.rho, 'x');

% gravity field of real model
[g_obs, fig_grav_obs] = calculate_gravity_field(Model_real.rho, rec_g);
% figname = [output_path,'/obs.gravityrecordings.png'];
% mtit(fig_grav_obs, 'gravity field of real model');
% print(fig_grav_obs, '-dpng', '-r400', figname);
% close(fig_grav_obs);

%% source time function
sEventInfoUnfilt = prepare_stf();
nsrc = length(sEventInfoUnfilt);

% make sources frequency dependent
for ifr = 1:nfr
    disp(['FREQUENCY NR. ',num2str(ifr),'/',num2str(nfr), ...
          '. fmin=',num2str(f_minlist(ifr)),', fmax=',num2str(f_maxlist(ifr))]);
    
    % get frequencies
    sObsPerFreq(ifr).f_max = f_maxlist(ifr);
    sObsPerFreq(ifr).f_min = f_minlist(ifr);
    
    % filter stf per src & per component
    sEventInfo = sEventInfoUnfilt;
    for isrc = 1:nsrc
        comps = fieldnames(sEventInfoUnfilt(isrc).stf);
        for icomp = 1:length(comps)
            stf = sEventInfoUnfilt(isrc).stf.(comps{icomp});
            t = sEventInfoUnfilt(isrc).t;
            stf = butterworth_lp(stf,t,3,sObsPerFreq(ifr).f_max,'silent');
            stf = butterworth_hp(stf,t,3,sObsPerFreq(ifr).f_min,'silent');
            sEventInfo(isrc).stf.(comps{icomp}) = stf; clearvars stf;
        end
    end
    
    % run forward per source
    [sEventObs, fig_seis] = run_forward_persource(Model_real, sEventInfo, 'yesplot');
    
    % plot & save each figure fig_seis
    for isrc = 1:nsrc
        titel = ['src ', num2str(isrc),' - observed seismograms freq range: ', ...
            num2str(sObsPerFreq(ifr).f_min), '-',num2str(sObsPerFreq(ifr).f_max), ' Hz'];
        mtit(fig_seis(isrc), titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = [output_path,'/obs.seis.fmax-',num2str(sObsPerFreq(ifr).f_max,'%.2e'),...
                               '.src-',num2str(isrc,'%02d'),'.png'];
        print(fig_seis(isrc),'-dpng','-r400',figname);
        close(fig_seis(isrc));
    end
    
    sObsPerFreq(ifr).sEventInfo = sEventInfo;
    sObsPerFreq(ifr).sEventObs = sEventObs;
end

t_obs = sObsPerFreq(ifr).sEventObs(1).t;

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = [output_path,'/obs.all-vars.mat'];
save(savename, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');

close all;

end

function Model_real = checkargs(args)

nargs = length(args);

if nargs ~= 1
    error('wrong input to prepare_obs !')
else
    if isnumeric(args{1})
        modelnr = args{1};
        Model_real = update_model(modelnr);
    elseif isstruct(args{1})
        Model_real = args{1};
    end
    
end


end