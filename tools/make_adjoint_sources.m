%==========================================================================
% compute and store adjoint sources -- one direction at a time! (x,y or z)
%
% input:
%-------
% v:            synthetic velocity seismograms
% v_0:          observed velocity seismograms
% t:            time axis
% output:       determines output type: 'displacement' or 'velocity'
% measurement:  'waveform_difference' for L2 waveform difference
%               'cc_time_shift' for cross-correlation time shift
% direction:    the appendix that gets added to each of the source time
%               function names. This is necessary for example if you want 
%               to calculate the source time functions for directions x, y,
%               z. The way it works in the code: 
%               for x '_1';   for y '_2';   for z '_3'
% mode:         'manual' or 'auto' (picking the seismogram)
%
% When v_0, i.e. the observed velocity seismograms, are set to zero, 
% the code performs data-independent measurements. 
%
% To get 'banana-doughnut kernels', i.e. travel time sensitivity kernels,
% use for measurement: 'cc_time_shift' and for output 'vel'.
%
% output:
%--------
% misfit:       The total misfit, calculated according to the current
%               chosen misfit functional.
% adjoint_stf:  The adjoint source time function for the specifict
%               direction. Time-reversed!!
%
% 'output' in terms of files written:
% adjoint sources are saved into ../input/sources/adjoint 
% 
%==========================================================================

function [adjoint_stf, misfit] = make_adjoint_sources(v,v_0,t,veldis, ...
                                                      measurement, ...
                                                      direction, mode)
%%
%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

path(path,'../input/');
path(path,'../code/propagation/');
path(path,'./misfits/')

% % delete previous instances of files
% delete(['../input/sources/adjoint/src*',direction])


input_parameters;
nrec = size(v,1);

% fid_loc=fopen([adjoint_source_path 'source_locations'],'w');

nt=length(t);
% length(v(1,:))
dt=t(2)-t(1);

misfit=0.0;

% initialise source time function output variable
adjoint_stf = zeros(nrec,nt); % adapt this to three dimensions when working!

%- convert to displacement if wanted ------------------------------------------

if strcmp(veldis,'displacement')
    
    %     u=zeros(nrec,nt);
    u=cumsum(v,2)*dt;
    v=u;
    u_0 = cumsum(v_0,2)*dt;
    v_0 = u_0;
elseif not(strcmp(veldis,'velocity'))
    error('ERRORRRR your output variable is not ''displacement'' or ''velocity''');
end


if strcmp(mode,'manual')
% reply = input('Do you want to manually pick the seismograms? ','s');
% if (strcmp(reply,'yes') || strcmp(reply,'y'))
    disp 'Manual labour... here we go!'
end
    
    
%==========================================================================
%% march through the various recordings -----------------------------------
%==========================================================================

% pick_data = figure;
adjoint_source = figure;

for n=1:nrec
    
    if strcmp(mode,'manual')
        fprintf(1,'station number %d\n',n)
    end
    
    %- plot traces --------------------------------------------------------
    
    figure(adjoint_source);
    clf;
    subplot(3,1,1)
    plot(t,v(n,:),'k')
    hold on
    plot(t,v_0(n,:),'r')
    plot(t,v(n,:)-v_0(n,:),'k--','LineWidth',2)
    hold off
    
    title(['receiver ' num2str(n) ' ,synth - black, obs - red, diff - dashed'])
    xlabel('t [s]')
    ylabel([veldis,' [m]'])
    
    
    
    %- select time windows and taper seismograms --------------------------
    if strcmp(mode,'auto')
%     if (strcmp(reply,'no') || strcmp(reply,'n'))
        left = t(1);
        right = t(end);
    elseif strcmp(mode,'manual')
%     elseif (strcmp(reply,'yes') || strcmp(reply,'y'))
        disp('select left window');
        [left,~]=ginput(1);
        disp('select right window');
        [right,~]=ginput(1);
    end
    
    width=t(end)/10;                                                    % why divide by ten??
    if (left < width) || (right > t(end) - width)
        left = min(left,width);
        right = max(right,t(end)-width);
    end
    v(n,:)=taper(v(n,:),t,left,right,width);
    v_0(n,:)=taper(v_0(n,:),t,left,right,width);
    
    %- compute misfit and adjoint source time function --------------------
    
    if strcmp(measurement,'waveform_difference')
        [misfit_n,adstf_nonreversed]=waveform_difference(v(n,:),v_0(n,:),t);
    elseif strcmp(measurement,'cc_time_shift')
        [misfit_n,adstf_nonreversed]=cc_time_shift(v(n,:),v_0(n,:),t);
    end
    
    misfit=misfit+misfit_n;
    
    
    %- plot adjoint source before time reversal -----------------------
    
    figure(adjoint_source);
    subplot(3,1,2);
    plot(t,adstf_nonreversed,'k')
    xlabel('t [s]')
    title(['adjoint source (', veldis, ' seismograms) before time reversal'])
    
    
    
    %- correction of adjoint source time function for velocity output -
    %  (because the signal is time-reversed, +velocity in the positive
    %  time direction becomes -velocity in the negative time direction)
    if strcmp(veldis,'velocity')
        adstf_nonreversed =-1 * adstf_nonreversed;
    end
    
    
    %- plot adjoint source after time reversal -----------------------
    figure(adjoint_source);
    subplot(3,1,3);
    plot(t,fliplr(adstf_nonreversed),'k')
    xlabel('t [s]')
    title(['adjoint source (', veldis, ' seismograms) after time reversal'])
    
    
    if strcmp(mode,'auto')
%     if (strcmp(reply,'no') || strcmp(reply,'n'))
        pause(0.05)
    elseif strcmp(mode,'manual')
%     elseif (strcmp(reply,'yes') || strcmp(reply,'yes'))
        pause(1.0)
    end
    
    %- save stf to adjoint_stf -- and time-reverse!! (using fliplr)
    
    adjoint_stf(n,:) = fliplr(adstf_nonreversed);
    
    
    
    
    %- write adjoint source locations to file -----------------------------
    
%     fprintf(fid_loc,'%g %g\n',rec_x(n),rec_z(n));
    
%     %- write source time functions ----------------------------------------
%     %  WITH time reversal!!!!!
    fn=[adjoint_source_path 'src_' num2str(n) direction];
    fid_src=fopen(fn,'w');
    for k=1:nt
        fprintf(fid_src,'%g\n',adstf_nonreversed(nt-k+1));
    end
    fclose(fid_src);
    
    
    
    %     end
end



% % save stf variable to a .mat file
% filename=['../input/sources/adjoint/adjoint_sources',direction];
% save(filename,'adjoint_stf','-v7.3')

%==========================================================================
%- clean up ---------------------------------------------------------------
%==========================================================================

% fclose(fid_loc);

end