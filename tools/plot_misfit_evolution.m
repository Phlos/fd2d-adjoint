function fig_misfit = plot_misfit_evolution(misfit_seis, misfit_g, misfit)

% plot the misfit evolution for different data together & separately.

set_figure_properties_bothmachines;


for i=1:length(misfit_seis)
    misf_seis(i) = misfit_seis(i).total;
    misf_grav(i) = misfit_g(i).total;
end

% set figure position
fig_misfit = figure;
set(fig_misfit,'OuterPosition',pos_misfit) 



%- plot various misfits:

% total
subplot(3,2,1);
plot(misfit,'k');
title('misfit evolution');
subplot(3,2,2);
semilogy(misfit,'k');
title('misfit evolution - semilogarithmic in misfit');

% seismic
subplot(3,2,3);
plot(misf_seis,'b');
title('misfit evolution (seis)');
subplot(3,2,4);
semilogy(misf_seis,'b');
title('misfit evolution (seis) - semilog in misfit');

% gravity
subplot(3,2,5);
plot(misf_grav,'r');
title('misfit evolution (grav)');
subplot(3,2,6);
semilogy(misf_grav,'r');
title('misfit evolution (grav) - semilog in misfit');

end