%- compute L2 waveform difference -----------------------------------------
%
% function [misfit,adsrc]=waveform_difference(u,u_0,t)
%
% u: synthetic displacement seismogram
% u_0: observed displacement seismogram
% t: time axis

function [misfit,adstf]=waveform_difference(u,u_0,t)

adstf=u-u_0;
misfit=sum(adstf.*adstf)*(t(2)-t(1));