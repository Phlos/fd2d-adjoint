%==========================================================================
% compute and store adjoint sources -- one direction at a time! (x,y or z)
%
% function misfit=make_adjoint_sources(u,u_0,t,veldis,measurement)
%
% input:
%-------
% v:            synthetic velocity seismograms
% v_0:          observed displacement seismograms
% t:            time axis
% veldis:       determines output type: 'displacement' or 'velocity'
% measurement:  'waveform_difference' for L2 waveform difference
%               'cc_time_shift' for cross-correlation time shift
% appendix:     the appendix that gets added to each of the source time
%               function names. This is necessary for example if you want 
%               to calculate the source time functions for directions x, y,
%               z. The way it works in the code: 
%               for x '_1';   for y '_2';   for z '_3'
%
% When v_0, i.e. the observed velocity seismograms, are set to zero, 
% the code performs data-independent measurements. This results in
% sensitivity kernels (i.e. the sensitivity of a certain observable -
% travel time, amplitude, ..., to a (infinitesimal) change in a model
% parameter.
%
% output:
%--------
% misfit:       The total misfit, calculated according to the current
%               chosen misfit functional.
%
% 'output' in terms of files written:
% adjoint sources are saved into ../input/sources/adjoint 
% 
%==========================================================================

function misfit=make_adjoint_sources(v,v_0,t,veldis,measurement,appendix)
%%
%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

path(path,'../input/');
path(path,'../code/propagation/');
path(path,'./misfits/')

input_parameters;
nrec = length(rec_x);

fid_loc=fopen([adjoint_source_path 'source_locations'],'w');

nt=length(t)
dt=t(2)-t(1)

misfit=0.0;

% initialise source time function save variable
adjoint_stf = zeros(nrec,nt); % adapt this to three dimensions when working!

%- convert to displacement if wanted ------------------------------------------

if strcmp(veldis,'displacement')
    
%     u=zeros(nrec,nt);
    u=cumsum(vel,2)*dt;
    v=u;
elseif not(strcmp(veldis,'velocity'))
    error('ERRORRRR your veldis variable is not ''displacement'' or ''velocity''');
end



%==========================================================================
%% march through the various recordings -----------------------------------
%==========================================================================

pick_data = figure;
adjoint_source = figure;

for n=1:nrec

        fprintf(1,'station number %d\n',n)
        
        %- plot traces --------------------------------------------------------
        
        figure(pick_data);
%         subplot(3,1,dir)
        plot(t,v(n,:),'k')
        hold on
        plot(t,v_0(n,:),'r')
        plot(t,v(n,:)-v_0(n,:),'k--')
        hold off
        
        title(['receiver ' num2str(n) ' ,synth - black, obs - red, diff - dashed'])
        xlabel('t [s]')
        ylabel('displacement [m]')
        
        %- select time windows and taper seismograms --------------------------
        
        disp('select left window');
        [left,~]=ginput(1);
        disp('select right window');
        [right,~]=ginput(1);
        
        width=t(end)/10;
        v(n,:)=taper(v(n,:),t,left,right,width);
        v_0(n,:)=taper(v_0(n,:),t,left,right,width);
        
        %- compute misfit and adjoint source time function --------------------
        
        if strcmp(measurement,'waveform_difference')
            [misfit_n,adstf]=waveform_difference(v(n,:),v_0(n,:),t);
        elseif strcmp(measurement,'cc_time_shift')
            [misfit_n,adstf]=cc_time_shift(v(n,:),v_0(n,:),t);
        end
        
        misfit=misfit+misfit_n;
        
        %- correct adjoint source time function for velocity measurement ------
        % ??????
%         if strcmp(veldis,'vel')
%             adstf(1:nt-1)=-diff(adstf)/dt;
%         end
        
        
        %- plot adjoint source before time reversal ---------------------------
        
        figure(adjoint_source);
        plot(t,adstf,'k')
        xlabel('t [s]')
        title(['adjoint source (', veldis, 'seismograms) before time reversal'])
        pause(1.0)
        
        %- write adjoint source locations to file -----------------------------
        
        fprintf(fid_loc,'%g %g\n',rec_x(n),rec_z(n));
        
        %- write source time functions ----------------------------------------
        %  WITH time reversal!!!!!
        fn=[adjoint_source_path 'src_' num2str(n) appendix];
        fid_src=fopen(fn,'w');
        for k=1:nt
            fprintf(fid_src,'%g\n',adstf(nt-k+1));
        end
        fclose(fid_src);
%     end
end

% save stf variable to a .mat file
save('../input/sources/adjoint/adjoint_sources','adjoint_stf','-v7.3')

%==========================================================================
%- clean up ---------------------------------------------------------------
%==========================================================================

fclose(fid_loc);