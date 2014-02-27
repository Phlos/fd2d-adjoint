%==========================================================================
% input parameters for interferometry
%==========================================================================

%- frequency sampling in Hz -----------------------------------------------
%- The sampling should be evenly spaced for the inverse F transform.-------
%- The sampling must also be sufficiently dense in order to avoid artefacts
%- on the positive time axis in the time-domain source function. This can 
%- be checked with "plot_correlation_source_function".

%- It is sufficient to consider the positive frequency axis. The
%- frequencies must start with 0.

f_sample=0.000:0.002:0.200;