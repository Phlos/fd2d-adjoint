function plot_vrec_to_adjointstf(t,v_rec,stf)

% to plot the original signal and the adjoint signal in one plot

figure;

% original recorded velocity seismogram
subplot(3,1,1);
plot(t,v_rec);
title('original seismograms');


% adjoint source, NON time reversed (= tapered velocity seismogram)
subplot(3,1,2);
plot(t,fliplr(stf));
title('non time-reversed');

% adjoint source, time reversed
subplot(3,1,3);
plot(t,stf);
title('time-reversed');

end