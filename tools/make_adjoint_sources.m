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
% appendix:     the appendix that gets added to each of the source time
%               function names. This is necessary for example if you want 
%               to calculate the source time functions for directions x, y,
%               z. The way it works in the code: 
%               for x '_1';   for y '_2';   for z '_3'
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
%
% 'output' in terms of files written:
% adjoint sources are saved into ../input/sources/adjoint 
% 
%==========================================================================

function [adjoint_stf, misfit] = make_adjoint_sources(v,v_0,t,output,measurement,direction)
%%
%==========================================================================
%- initialisations --------------------------------------------------------
%==========================================================================

path(path,'../input/');
path(path,'../code/propagation/');
path(path,'./misfits/')

% delete previous instances of files
 delete(['../input/sources/adjoint/src*',direction])

input_parameters;
nrec = length(rec_x);

fid_loc=fopen([adjoint_source_path 'source_locations'],'w');

nt=length(t);
dt=t(2)-t(1);

misfit=0.0;

% initialise source time function save variable
adjoint_stf = zeros(nrec,nt); % adapt this to three dimensions when working!

%- convert to displacement if wanted ------------------------------------------

if strcmp(output,'displacement')
    
%     u=zeros(nrec,nt);
    u=cumsum(vel,2)*dt;
    v=u;
elseif not(strcmp(output,'velocity'))
    error('ERRORRRR your output variable is not ''displacement'' or ''velocity''');
end



%==========================================================================
%% march through the various recordings -----------------------------------
%==========================================================================

% pick_data = figure;
adjoint_source = figure;

for n=1:nrec

        fprintf(1,'station number %d\n',n)
        
        %- plot traces --------------------------------------------------------
        
        figure(adjoint_source);
        subplot(2,1,1)
        plot(t,v(n,:),'k')
        hold on
        plot(t,v_0(n,:),'r')
        plot(t,v(n,:)-v_0(n,:),'k--')
        hold off
        
        title(['receiver ' num2str(n) ' ,synth - black, obs - red, diff - dashed'])
        xlabel('t [s]')
        ylabel('velocity [m]')
        
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
%         if strcmp(veldis,'displacement')
%             adstf(1:nt-1)=-diff(adstf)/dt;
%         end
        
        
        %- plot adjoint source before time reversal ---------------------------
        
        figure(adjoint_source);
        subplot(2,1,2);
        plot(t,adstf,'k')
        xlabel('t [s]')
        title(['adjoint source (', output, 'seismograms) before time reversal'])
        pause(1.0)
        
        %- write adjoint source locations to file -----------------------------
        
        fprintf(fid_loc,'%g %g\n',rec_x(n),rec_z(n));
        
        %- write source time functions ----------------------------------------
        %  WITH time reversal!!!!!
        fn=[adjoint_source_path 'src_' num2str(n) direction];
        fid_src=fopen(fn,'w');
        for k=1:nt
            fprintf(fid_src,'%g\n',adstf(nt-k+1));
        end
        fclose(fid_src);
        
        %- save source time functions to adjoint_stf
        adjoint_stf(n,:) = adstf;
        
%     end
end

% save stf variable to a .mat file
filename=['../input/sources/adjoint/adjoint_sources',direction];
save(filename,'adjoint_stf','-v7.3')

%==========================================================================
%- clean up ---------------------------------------------------------------
%==========================================================================

fclose(fid_loc);