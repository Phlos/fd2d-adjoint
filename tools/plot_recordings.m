% function plot_recordings(u,t,mode)
%
% v:        velocity recordings (dimensions (n_receivers, nt) )
% t:        time axis
% veldis:   the mode of the seismograms that we plot: 'displacement' or
% 'velocity'

function plot_recordings(vel,t,veldis)

%==========================================================================
%- plot recordings, ordered according to distance from the first source ---
%==========================================================================

%- initialisations and parameters -----------------------------------------

spacing=1.5;
sort=0;

%- read input -------------------------------------------------------------

path(path,'../input/');
input_parameters;

%- make distance vector and sort ------------------------------------------

if (sort==1)

    d=sqrt((rec_x-src_x).^2+(rec_z-src_z).^2);
    [dummy,idx]=sort(d);

else
    
    idx=1:length(rec_x);
    
end

%- convert to displacement if wanted ------------------------------------------

if strcmp(veldis,'displacement')
    nt=length(t);
    u=zeros(length(rec_x),nt);
    
    u = cumsum(vel,2)*dt;
    
    % now put the displacement values back in the vel variable so that it
    % can be plotted.
    vel=u;
elseif (not(strcmp(veldis,'velocity')))
    error('ERRORR your veldis input variable is wrong. Eejit.');
end

%- plot recordings with ascending distance from the first source ----------

recordings = figure;
set(gca,'FontSize',20)
hold on

for k=1:length(rec_x)
    
    m=max(abs(vel(idx(k),:)));
    % traces om en om blauw en zwart kleuren
    if mod(k,2)
        plot(t,spacing*(k-1)+vel(idx(k),:)/m,'k','LineWidth',1)
    else
        plot(t,spacing*(k-1)+vel(idx(k),:)/m,'b','LineWidth',1)
    end
    
    % afstanden geven in m of km
    if (max(rec_x)<=1000)
        text(min(t)-(t(end)-t(1))/6,spacing*(k-1)+0.3,['x=' num2str(rec_x(idx(k))) ' m, z=' num2str(rec_z(idx(k))) ' m'],'FontSize',14)
    else
        text(min(t)-(t(end)-t(1))/6,spacing*(k-1)+0.3,['x=' num2str(rec_x(idx(k))/1000) ' km, z=' num2str(rec_z(idx(k))/1000) ' km'],'FontSize',14)
    end
    
end

% text on figure: velocity or displacement seismogram
figure(recordings);
title([veldis, ' seismograms'])

xlabel('time [s]','FontSize',20);
ylabel('normalised traces','FontSize',20);
axis([min(t)-(t(end)-t(1))/5 max(t)+(t(end)-t(1))/10 -1.5 spacing*length(rec_x)+1])