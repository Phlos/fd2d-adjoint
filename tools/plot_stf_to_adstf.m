function fig_stftoadstf = plot_stf_to_adstf(stf, vel, adstf, t)

% This function plots a figure that shows how the source time function
% changes into the adjoint source time function
%
% SYNTAX:
% fig_stftoadstf = plot_stf_to_adstf(stf, disp, vel, adstf, t);
%
% INPUT:
% - stf:    the original source time function
% - disp:   the displacement at the station
% - vel:    the velocity at the station
% - adstf:  the final adjoint source time function at the station
% - t:      the time axis
%
% OUTPUT:
% - figure  plotting stf, disp, vel, adstf

%- preparation
dt = t(2)-t(1);
disp = cumsum(vel)*dt;

%- actual plotting

fig_stftoadstf = figure;

subplot(2,2,1)
plot(t,stf);
xlabel('time (s)');
ylabel({'stf amplitude';' unknown units'});
title('source-time function');

subplot(2,2,2)
plot(t,disp);
xlabel('time (s)');
ylabel('displacement (m)');
title('displacement seismogram');

subplot(2,2,3)
plot(t,vel);
xlabel('time (s)');
ylabel('velocity  (m/s)');
title('velocity seismogram');

subplot(2,2,4)
plot(t,adstf);
xlabel('time (s)');
ylabel('adjoint stf, unknown units');
title('adjoint source-time function');

end