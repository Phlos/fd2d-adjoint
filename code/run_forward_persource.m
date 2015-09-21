function [sEventRec, fig_seis] = run_forward_persource(Model, sEventInfo, varargin)

% wrapper to run forward for each source consecutively



%% prep
input_parameters;
% project folder
output_path = ['./output/',project_name,'/'];
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
            if ~(exist([output_path,'/fwd_temp/'], 'dir')); 
                mkdir([output_path,'/fwd_temp/']); 
            end
            prevmsg = sprintf( 'saving u_fw, v_fw output to file...');
            fprintf(prevmsg);
            varInfo = whos('u_fw', 'v_fw');
            size_together = varInfo(1).bytes + varInfo(2).bytes;
            if size_together < 2.1472e9
                save([output_path,'/fwd_temp/forwardfield.src-',num2str(isrc),'.mat'], 'u_fw', 'v_fw', '-v6');
            else
                warning('MAT-file size exceeds maximum for -v6 -- will save as v7.3');
                save([output_path,'/fwd_temp/forwardfield.src-',num2str(isrc),'.mat'], 'u_fw', 'v_fw', '-v7.3');
            end
            reverseStr = repmat(sprintf('\b'), 1, length(prevmsg));
            fprintf(reverseStr);
            clearvars u_fw v_fw;
        end

        if strcmp(plotornot,'yesplot')
            % determine how many seismograms are actually plotted
            recs=nrec;
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

% determine some output:
% - whether seisdif figure can use obs
% - whether the forward wavefield should be saved
% - whether plots should be made for seisdif etc.

% default values:
obsPresent = 'no';
sEventObs = NaN;
saveFwdFields = 'nosavefields';
plotornot = 'noplot';

for ii = 1:numel(args)
    if isstruct(args{ii}) && isfield(args{ii}, 'vel')
        obsPresent = 'yes';
        sEventObs = args{1};
    elseif ischar(args{ii})
        if any(strcmp(args{ii}, {'yesplot'; 'noplot'}))
            plotornot = args{ii};
        elseif any (strcmp(args{ii}, {'yessavefields'; 'nosavefields'}))
            saveFwdFields = args{ii};
        end
    else
        error(['Obs input nr',num2str(ii), ': ',args{ii}, ' to run_foward_persource is ambiguous']);
    end
end
    
end