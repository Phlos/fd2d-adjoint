function [sEventRec, fig_seis] = run_forward_persource(Model, sEventInfo, varargin)

% wrapper to run forward for each source consecutively

%% prep
nsrc = length(sEventInfo);

% determine if v_obs is supplied (in order to plot seisdif)
[obsPresent, sEventObs, saveFwdFields, plotornot] = checkargs(varargin(:));

% prepare zeros sources
for isrc = 1:nsrc
    stf_zero{isrc} = make_seismogram_zeros(sEventInfo(isrc).stf);
end

%% actual loop over sources

    % run fwd per source (other sources = zeros)
    for isrc = 1:nsrc
        disp(['Running forward wave propagation - src. nr. ',num2str(isrc),'/',num2str(nsrc)]);
        
        % add a single nonzero source for isrc
        stf = stf_zero;
        stf{isrc} = sEventInfo(isrc).stf;
        
        % run actual fwd wave propagation per src per freq
        [vel,t,u_fw,v_fw,rec_x,rec_z] = run_forward(Model, stf);

        % save forward fields to file
        if strcmp(saveFwdFields, 'yessavefields')
            disp 'saving u_fw, v_fw output to file...'
            save(['../output/forwardfield.src-',num2str(isrc),'.mat'], 'u_fw', 'v_fw', '-v6');
            clearvars u_fw v_fw;
        end

        if strcmp(plotornot,'yesplot')
            % determine how many seismograms are actually plotted
            if length(vel) > 8; recs = [2:2:length(vel)]; end
            
            if strcmp(obsPresent, 'yes');
                vobs = sEventObs(isrc).vel;
                fig_seis(isrc) = plot_seismogram_difference(vel, vobs, t, recs);
            else
                vobs = make_seismogram_zeros(vel);
                fig_seis(isrc) = plot_seismogram_difference(vel, vobs, t, recs, 'nodiff');
            end
        else
            fig_seis = NaN;
        end
        
        
    
        % write obs into freq & source gather variable
        sEventRec(isrc).vel = vel;
        sEventRec(isrc).t = t;
        sEventRec(isrc).rec_x = rec_x;
        sEventRec(isrc).rec_z = rec_z;
        
    end

end

function [obsPresent, sEventObs, saveFwdFields, plotornot] = checkargs(args)

% determine whether seisdif figure can use obs

% default values:
obsPresent = 'no';
sEventObs = NaN;
saveFwdFields = 'nosavefields';
plotornot = 'noplot';

% narg = length(args);
% 
% if narg == 0
%     obsPresent = 'no';
%     sEventObs = NaN;
% elseif narg == 1 
%     if isfield(args{1}, 'vel')
%         obsPresent = 'yes';
%         sEventObs = args{1};
%     elseif ischar(args{1})
%         saveFwdFields = args{1};
%     end
% elseif narg == 2 && isfield(args{1}, 'vel') && ischar(args{2})
%     obsPresent = 'yes';
%     sEventObs = args{1};
%     saveFwdFields = args{2};
% else
%     error('Obs input to run_foward_persource is ambiguous');
% end

for ii = 1:numel(args)
    if isstruct(args{ii}) && isfield(args{ii}, 'vel')
        obsPresent = 'yes';
        sEventObs = args{1};
    elseif ischar(args{ii})
        if any(strcmp(args{ii}, {'yesplot'; 'noplot'}))
            plotornot = args{ii};
        elseif any (strcmp(args{ii}, {'yessavefields'; 'nosavefields'}))
            saveFwdfields = args{ii};
        end
    else
        error(['Obs input ',args{ii}, ' to run_foward_persource is ambiguous']);
    end
end
    
end