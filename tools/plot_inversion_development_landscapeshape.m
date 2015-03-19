function fig_inv2 = plot_inversion_development_landscapeshape(InvProps, imax)

% plot misfit development with nonlinearity
% SYNTAX:
% fig_inv2 = plot_inversion_development_landscapeshape(InvProps, imax);

%% prepare

set_figure_properties_bothmachines;

% variables
misfit = InvProps.misfit;
misfit_seis = InvProps.misfit_seis;
misfit_g = InvProps.misfit_g;
step = InvProps.step;

niter = length(misfit);
iters = 2:niter;

for i = 1:niter
misfitseis(i) = misfit_seis{i}.normd;
end
for i = 1:niter
misfitgrav(i) = misfit_g(i).normd;
end

stepcumsum = cumsum(abs(step(1:imax-1)));

% plot prep
fig_inv2 = figure;
pos_invplot2 = [170         232         829        1411];
set(fig_inv2, 'OuterPosition', pos_invplot2)

%% plot
% misfit development, semilogy
subplot(3,1,1)
h = semilogy(iters,misfit(2:imax), '-k', iters,misfitseis(2:imax), '-r', iters, misfitgrav(2:imax), '-b');
set(h(1), 'LineWidth', 2)
set(h(2), 'LineWidth', 1)
set(h(3), 'LineWidth',1)
legend('total', 'seis', 'grav')
title('misfit development (norm''d)');
xlabel('iteration no.');
ylabel('misfit');

% misfit versus step - an indication of how irregular the landscape is
subplot(3,1,2)
% h2 = plot(stepcumsum / stepcumsum(end) * 100,misfit(2:100), 'k-x');
h2 = semilogy(stepcumsum,misfit(2:imax), 'k-x', ...
    stepcumsum, misfitseis(2:imax), 'r-x', ...
    stepcumsum, misfitgrav(2:imax), '-bx');
xlim([min([stepcumsum(1:imax-1), 0]) max(stepcumsum(1:imax-1))]);
set(h2(1), 'LineWidth', 2)
set(h2(2), 'LineWidth', 1)
set(h2(3), 'LineWidth',1)
title('misfit vs. step length');
xlabel('distance along path in misfit landscape');
ylabel('misfit');

% step length -- just checking if anything is correlated
subplot(3,1,3)
h3 = semilogy(stepcumsum,step(1:imax-1), 'k');
xlim([0 max(stepcumsum)]);
set(h3, 'LineWidth', 1)
title('step length');
xlabel('distance along path in misfit landscape');
ylabel('step length');

end