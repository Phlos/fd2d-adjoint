%==========================================================================
% compute and store adjoint sources -- one direction at a time! (x,y or z)
%
% function misfit=make_adjoint_sources(u,u_0,t,veldis,measurement)
%
% input:
%-------
% u:            synthetic displacement seismograms
% u_0:          observed displacement seismograms
% t:            time axis
% veldis:       'dis' for displacements, 'vel' for velocities
% measurement:  'waveform_difference' for L2 waveform difference
%               'cc_time_shift' for cross-correlation time shift
%
% When u_0, i.e. the observed displacement seismograms, are set to zero, 
% the code performs data-independent measurements. 
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

function misfit=make_adjoint_sources(u,u_0,t,veldis,measurement)
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

%- convert to velocity if wanted ------------------------------------------

if strcmp(veldis,'vel')
    
    v=zeros(nrec,nt);
    
    for k=1:nrec
        v(k,1:nt-1)=diff(u(k,:))/(t(2)-t(1));
        v(k,nt)=0.0;
    end
   
    u=v;
    
end

%%
%==========================================================================
%- march through the various recodings ------------------------------------
%==========================================================================

pick_data = figure;
adjoint_source = figure;

for n=1:nrec
%     for dir = 1:3
%         if (dir==1)
%             u=ux;
%             u_0=ux_0;
%         elseif (dir==2)
%             u=uy;
%             u_0=uy_0;
%         elseif (dir==3)
%             u=uz;
%             u_0=uz_0;
%         else
%             disp 'ERROR direction exceeds dimensions'
%         end
        fprintf(1,'station number %d\n',n)
        
        %- plot traces --------------------------------------------------------
        
        figure(pick_data);
%         subplot(3,1,dir)
        plot(t,u(n,:),'k')
        hold on
        plot(t,u_0(n,:),'r')
        plot(t,u(n,:)-u_0(n,:),'k--')
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
        u(n,:)=taper(u(n,:),t,left,right,width);
        u_0(n,:)=taper(u_0(n,:),t,left,right,width);
        
        %- compute misfit and adjoint source time function --------------------
        
        if strcmp(measurement,'waveform_difference')
            [misfit_n,adstf]=waveform_difference(u(n,:),u_0(n,:),t);
        elseif strcmp(measurement,'cc_time_shift')
            [misfit_n,adstf]=cc_time_shift(u(n,:),u_0(n,:),t);
        end
        
        misfit=misfit+misfit_n;
        
        %- correct adjoint source time function for velocity measurement ------
        
        if strcmp(veldis,'vel')
            adstf(1:nt-1)=-diff(adstf)/dt;
        end
        
        
        %- plot adjoint source before time reversal ---------------------------
        
        figure(adjoint_source);
        plot(t,adstf,'k')
        xlabel('t [s]')
        title('adjoint source before time reversal')
        pause(1.0)
        
        %- write adjoint source locations to file -----------------------------
        
        fprintf(fid_loc,'%g %g\n',rec_x(n),rec_z(n));
        
        %- write source time functions ----------------------------------------
        %  WITH time reversal!!!!!
        fn=[adjoint_source_path 'src_' num2str(n)];
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