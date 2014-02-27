%- compute cross-correlation time shift -----------------------------------
%
% function [misfit,adsrc]=cc_time_shift(u,u_0,t)
%
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% t: time axis

function [misfit,adstf]=cc_time_shift(u,u_0,t)

%- initialisations --------------------------------------------------------

path(path,'tf/');

%- compute time shift -----------------------------------------------------

if sum(u_0==0)==length(t)
    T=1.0;
else
    [cc,t_cc]=cross_correlation(u_0,u,t);
    T=t_cc(find(max(cc)==cc));
end

%- compute adjoint source time function -----------------------------------

dt=t(2)-t(1);
nt=length(t);
v=zeros(1,nt);

v(1:nt-1)=diff(u)/dt;

adstf=v/(sum(v.*v)*dt);

%- compute misfit ---------------------------------------------------------

misfit=T^2/2.0;