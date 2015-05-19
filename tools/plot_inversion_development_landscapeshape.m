function fig_inv2 = plot_inversion_development_landscapeshape(InvProps, imax)

% plot misfit development with nonlinearity
% SYNTAX:
% fig_inv2 = plot_inversion_development_landscapeshape(InvProps, imax);

%% prepare

set_figure_properties_bothmachines;

% variables
misfit = InvProps.misfit;
misfitseis = InvProps.misfitseis;
misfitgrav = InvProps.misfitgrav;
step = InvProps.step;

niter = length(misfit);
iters = 1:niter;


stepcumsum = [0, cumsum(abs(step(1:imax-1)))];

% plot prep
fig_inv2 = figure;
pos_invplot2 = [170         232         829        1411];
set(fig_inv2, 'OuterPosition', pos_invplot2)

%% plot
% misfit development, semilogy
subplot(3,1,1)
% imax
% length(iters)
h = semilogy(iters(1:imax),misfit(1:imax), '-k', ...
             iters(1:imax),misfitseis(1:imax), '-r', ...
             iters(1:imax), misfitgrav(1:imax), '-b');
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
h2 = semilogy(stepcumsum(1:imax),misfit(1:imax), 'k-x', ...
    stepcumsum(1:imax), misfitseis(1:imax), 'r-x', ...
    stepcumsum(1:imax), misfitgrav(1:imax), '-bx');
% min([stepcumsum(1:imax), 0])
% max(stepcumsum(1:imax))
xlim([min([stepcumsum(1:imax), 0]) max(stepcumsum(1:imax))]);
set(h2(1), 'LineWidth', 2)
set(h2(2), 'LineWidth', 1)
set(h2(3), 'LineWidth',1)
title('misfit vs. step length');
xlabel('distance along path in misfit landscape');
ylabel('misfit');

% step length -- just checking if anything is correlated
subplot(3,1,3)
h3 = semilogy(stepcumsum(2:imax),step(1:imax-1), 'k');
xlim([0 max(stepcumsum)]);
set(h3, 'LineWidth', 1)
title('step length');
xlabel('distance along path in misfit landscape');
ylabel('step length');

end