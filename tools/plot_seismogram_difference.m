function fig_seisdif = plot_seismogram_difference(v_obs, v_rec, t, varargin)

% This function is (just like @ make_adjoint_sources) to plot the observed
% and recorded seismograms, and the difference between them.
%
% SYNTAX:
% fig_seisdif = plot_seismogram_difference(v_obs, v_rec, t)
%           and plot_seismogram_difference(v_obs, v_rec, t, 'yesdiff')
%               [plots seismograms from all stations obs, rec and diff]
% fig_seisdif = plot_seismogram_difference(v_obs, v_rec, t, 'nodiff')
%               [plots seismograms from all stations obs, rec - NO DIFF]
% fig_seisdif = plot_seismogram_difference(v_obs, v_rec, t, [recs])
%           and plot_seismogram_difference(v_obs, v_rec, t, [recs], 'yesdiff')
%               [plots seismograms obs, rec and diff for stations defined in [recs] ]
% fig_seisdif = plot_seismogram_difference(v_obs, v_rec, t, [recs], 'nodiff')
%               [plots seismograms obs, rec for stations defined in [recs] - NO DIFF]
%
% INPUT:
% - v_obs:  struct containing x and/or y and/or z traces of seismograms
% - v_rec:  struct containing x and/or y and/or z traces of seismograms. At
%           the very least, the same components as are present in v_obs
%           must be present in this struct, else errors will ensue.
% - t:      time axis.
%
% OUTPUT:
% - figure  plotting both sets of seismograms, plus the difference traces of
%   (v_rec - v_obs).
% - seisdif: figure handle of this figure

fig_seisdif = figure;
set_figure_properties_bothmachines;
set(fig_seisdif, 'OuterPosition', pos_seis);

[recs_given, recs_in, plot_diff] = check_args(varargin(:));


% number of components (x,y,z) for which seismograms have been recorded
ncomp = size(fieldnames(v_obs{1}), 1);

switch recs_given
    case 'yes'
        nrec = length(recs_in);
        recs = recs_in;
    case 'no'
        % number of receivers for which we have seismograms
        nrec = length(v_obs);
        recs = 1:nrec;
    otherwise
        error('wrong recs');
end

iplot = 0;
maks = 0;
for irec = recs
    comp = fieldnames(v_obs{irec});
    for icomp = 1:length(comp);
        
%         subplot(ncomp,1,icomp);
        iplot = iplot+1;
        assen(iplot) = subplot(nrec,ncomp,(irec-1)*ncomp + icomp);
        hold on
        plot(t,v_rec{irec}.(comp{icomp}),'k', ...
             t,v_obs{irec}.(comp{icomp}),'r--');
         if strcmp(plot_diff, 'yesdiff')
             plot(t,v_rec{irec}.(comp{icomp}) - v_obs{irec}.(comp{icomp}), 'b');
         end
        maks = max(maks, max(assen(iplot).YLim));
        %     plot(t,v_rec.(comp{1}) - v_obs.(comp{1}), 'b', 'LineWidth',2)
%         hold off
        
        if irec==1
            if strcmp(plot_diff, 'yesdiff')
                title({[comp{icomp},' component:']; 'synth - black, obs - red, (synth-obs) - blue'})
            else
                title({[comp{icomp},' component:']; 'synth - black, obs - red'})
            end
        end
        if irec==nrec
        xlabel('t [s]')
        end
        ylabel('v [m/s]')
        
%         icomp=icomp+1;
        
    end
end

linkaxes(assen, 'y');
assen(iplot).YLim = [-maks maks];

end

function [recs_given, recs, plot_diff] = check_args(args)

nargs = length(args);

switch nargs
    case 0
        plot_diff = 'yesdiff';
        recs_given = 'no';
        recs = NaN;
    case 1
        if ischar(args{1})
            plot_diff = args{1};
            recs_given = 'no';
            recs = NaN;
        else
            plot_diff = 'yesdiff';
            recs_given = 'yes';
            recs = args{1};
        end
    case 2
        recs_given = 'yes';
        recs = args{1};
        plot_diff = args{2};
    otherwise
        error('unknown number of input args')
end

end