function [Model_real, freqs, t_obs, props_obs, g_obs] = prepare_obs(output_path, varargin)

input_parameters;
nfr = length(f_maxlist);

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
[sources, t] = prepare_stf();
nsrc = length(sources);

% make sources frequency dependent
for ifr = 1:nfr
    
    % get frequencies
    freqs(ifr).f_max = f_maxlist(ifr);
    freqs(ifr).f_min = f_minlist(ifr);
    
    % filter stf per src & per component
    for isrc = 1:length(sources)
        comps = fieldnames(sources(isrc).stf);
        for icomp = 1:length(comps)
            stf = sources(isrc).stf.(comps{icomp});
            stf = butterworth_lp(stf,t,3,freqs(ifr).f_max,'silent');
            stf = butterworth_hp(stf,t,3,freqs(ifr).f_min,'silent');
            freqs(ifr).source(isrc).stf.(comps{icomp}) = stf; clearvars stf;
        end
    end
    
    % prepare zeros sources
    for isrc = 1:nsrc
        stf_zero{isrc} = make_seismogram_zeros(freqs(ifr).source(isrc).stf);
    end
    
    % run fwd per source (other sources = zeros)
    for isrc = 1:nsrc
        disp(['Making obs - freq nr. ',num2str(ifr),'/',num2str(nfr),' - src. nr ',num2str(isrc),'/',num2str(nsrc)]);
        
        % add a single nonzero source for isrc
        stf = stf_zero;
        stf{isrc} = freqs(ifr).source(isrc).stf;
        
        % run actual fwd wave propagation per src per freq
        [vobs,t_obs,~,~,~,~] = run_forward(Model_real, stf);

        v_0 = make_seismogram_zeros(vobs);
        fig_seis = plot_seismogram_difference(vobs, v_0, t_obs, 'nodiff');
        titel = [project_name,': observed seismograms freq range: ', ...
                 num2str(freqs(ifr).f_min), '-',num2str(freqs(ifr).f_max), ' Hz'];
        mtit(fig_seis, titel, 'xoff', 0.001, 'yoff', -0.05);
        figname = [output_path,'/obs.seis.fmax-',num2str(freqs(ifr).f_max,'%.2e'),'.png'];
        print(fig_seis,'-dpng','-r400',figname);
%         pause(1);
        close(fig_seis);
    
        % write obs into freq & source gather variable
        freqs(ifr).source(isrc).v_obs = vobs;
    end
end

% saving the obs variables
disp 'saving obs variables to matfile...'
savename = [output_path,'/obs.all-vars.mat'];
save(savename, 'freqs', 't_obs', 'Model_real', 'props_obs', 'g_obs', '-v7.3');

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