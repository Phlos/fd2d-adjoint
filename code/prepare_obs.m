function [Model_real, sObsPerFreq, t_obs, props_obs, g_obs] = prepare_obs(output_path, varargin)

input_parameters;
nfr = length(f_maxlist);


%% Model
[Model_real, onlyseiss] = checkargs(varargin(:));
if strcmp(onlyseiss,'yesonlyseis')
    onlyseis = true;
else
    onlyseis = false;
end

% Model_real = update_model(modelnr);

% % plotting the real model
% fig_mod_real = plot_model(Model_real, 'rhovsvp');
% mtit(fig_mod_real, 'real model -- rho-vs-vp parametrisation');
% figname = [output_path,'/obs.model.rhovsvp.png'];
% print(fig_mod_real, '-dpng', '-r400', figname);
% close(fig_mod_real);

% real - background model
if ~onlyseis
    fig_mod_diff = plot_model_diff(Model_real, bg_model_type, 'rhovsvp');
    mtit(fig_mod_diff, 'real - bg model -- rho-vs-vp parametrisation');
    figname = [output_path,'/obs.model-real-diff-bg.rhovsvp.png'];
    print(fig_mod_diff, '-dpng', '-r400', figname);
    close(fig_mod_diff);
end

% h.c. model properties for real model
props_obs = calculate_model_properties(Model_real.rho, 'x');

% gravity field of real model
    [g_obs, fig_grav_obs] = calculate_gravity_field(Model_real.rho, rec_g);
if ~onlyseis
    figname = [output_path,'/obs.gravityrecordings.png'];
    mtit(fig_grav_obs, 'gravity field of real model');
    print(fig_grav_obs, '-dpng', '-r400', figname);
    close(fig_grav_obs);
end

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
if ~onlyseis
    disp 'saving obs variables to matfile...'
    savename = [output_path,'/obs.all-vars.mat'];
    save(savename, 'sObsPerFreq', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v6');
end

close all;

end

function [Model_real, onlyseis] = checkargs(args)

nargs = length(args);

onlyseis = 'noonlyseis';

% if nargs ~= 1
%     error('wrong input to prepare_obs !')
% else
for ii = 1:nargs
    if isnumeric(args{ii})
        modelnr = args{ii};
        Model_real = update_model(modelnr);
    elseif isstruct(args{ii})
        Model_real = args{ii};
    elseif ischar(args{ii})
        onlyseis = args{ii};
    else
        error('wrong input to prepare_obs !')
    end
end
% end


end