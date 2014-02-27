% function plot_recordings(u,t,mode)
%
% u: displacement recordings
% t: time axis
% mode: 'dis' for displacement, 'vel' for velocity

function plot_recordings(u,t,veldis)

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

%- convert to velocity if wanted ------------------------------------------

if strcmp(veldis,'vel')
    nt=length(t);
    v=zeros(length(rec_x),nt-1);
    
    for k=1:length(rec_x)
        v(k,:)=diff(u(k,:))/(t(2)-t(1));
    end
    
    t=t(1:nt-1);
    u=v;
    
end

%- plot recordings with ascending distance from the first source ----------

recordings = figure;
set(gca,'FontSize',20)
hold on

for k=1:length(rec_x)
    
    m=max(abs(u(idx(k),:)));
    if mod(k,2)
        plot(t,spacing*(k-1)+u(idx(k),:)/m,'k','LineWidth',1)
    else
        plot(t,spacing*(k-1)+u(idx(k),:)/m,'b','LineWidth',1)
    end
    
    if (max(rec_x)<=1000)
        text(min(t)-(t(end)-t(1))/6,spacing*(k-1)+0.3,['x=' num2str(rec_x(idx(k))) ' m, z=' num2str(rec_z(idx(k))) ' m'],'FontSize',14)
    else
        text(min(t)-(t(end)-t(1))/6,spacing*(k-1)+0.3,['x=' num2str(rec_x(idx(k))/1000) ' km, z=' num2str(rec_z(idx(k))/1000) ' km'],'FontSize',14)
    end
    
end

xlabel('time [s]','FontSize',20);
ylabel('normalised traces','FontSize',20);
axis([min(t)-(t(end)-t(1))/5 max(t)+(t(end)-t(1))/10 -1.5 spacing*length(rec_x)+1])