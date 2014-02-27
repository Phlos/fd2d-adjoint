function s=ricker_wavelet(stdvar,t0,t)

%==========================================================================
% ricker wavelet = time derivative of a Gaussian as a function of time t
% with standard deviation sigma and time shift t0
%==========================================================================

s=-(t-t0)/stdvar^2;
s=s.*exp(-0.5*(t-t0).^2/stdvar^2);
