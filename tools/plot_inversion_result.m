function fig_invres = plot_inversion_result(InvProps, imax)

% plotting the results of the inversion in a good format.
% 3 x 2 subplots with:
% (left column:)
% - misfit(s) vs. iter
% - modeldifnorm vs. iter
% - step size vs. iter
% (right column:)
% - misfit vs. cumulative steplength
% - angle vs. misfit path? or iter?!

%% prepare

input_parameters;
set_figure_properties_bothmachines;

% variables
misfit = InvProps.misfit;
misfit_seis = InvProps.misfit_seis;
misfit_g = InvProps.misfit_g;
modeldifn = InvProps.modeldifn;
step = InvProps.step;
angletot = InvProps.angle.Ktotal;
angleseis = InvProps.angle.Kseis;
if strcmp(use_grav,'yes')
anglegrav = InvProps.angle.Kg;
else
    anglegrav = NaN*zeros(size(angleseis));
end

% some adaptations
niter = length(misfit); % the actual total number of iters run
imax = min(imax,niter);
iters = 1:niter;

for ii = 1:imax
    misfitseis(ii) = misfit_seis{ii}.normd;
    misfitgrav(ii) = misfit_g(ii).normd;
end

stepcumsum = [0, cumsum(abs(step(1:imax-1)))];

% % check output length
% disp(['length misfit: ', num2str(length(misfit))]);
% disp(['length seis misfit: ', num2str(length(misfitseis))]);
% disp(['length grav misfit: ', num2str(length(misfitgrav))]);
% disp(['length modeldifn: ', num2str(length(modeldifn))]);
% disp(['length step: ', num2str(length(step))]);
% disp(['length angletot: ', num2str(length(angletot))]);
% disp(['length angleseis: ', num2str(length(angleseis))]);
% disp(['length anglegrav: ', num2str(length(anglegrav))]);
% disp(['length iters: ', num2str(length(iters))]);
% disp(['length stepcumsum: ', num2str(length(stepcumsum))]);


%- prepare figure
fig_invres = figure;
set(fig_invres,'OuterPosition',pos_invres);


% s(1) = 

%% - left column: vs iter

% misfit development
subplot(9,2,[1,5])
Lmisfit = semilogy(iters(1:imax), misfit(1:imax), 'k',... 
        iters(1:imax), misfitseis(1:imax), '-r', ...
        iters(1:imax), misfitgrav(1:imax), '-b');
    set(Lmisfit(1), 'LineWidth', 2)
    set(Lmisfit(2), 'LineWidth', 1)
    set(Lmisfit(3), 'LineWidth',1)
%     h(1).LineWidth = 2;
% xlabel('iteration no.');
% ylabel('misfit (normalised)');
xlim([min(iters), iters(imax)]);
grid on
title('misfit (normalised)');
legend('all data combined', 'seismic data only', 'gravity data only', 'Location', 'Northeast');

% norm of |model - prev model|
% subplot(4,4,13)
subplot(9,2,[7 9])
% plot(modeldifn,'k');
if imax > 1
Lmdifn = semilogy(iters(2:imax), modeldifn(2:imax),'k');
end
xlim([min(iters), iters(imax)]);
grid on
set(Lmdifn, 'LineWidth', 1)
title('|current - previous model| / |1st model|');

% angle between consecutive kernels - vs. iter
subplot(9,2,[11 13])
if imax > 2
Langle = plot(iters(2:imax-1), angletot(2:imax-1) ./ pi .* 180, 'k',... 
        iters(2:imax-1), angleseis(2:imax-1) ./ pi .* 180, '-r', ...
        iters(2:imax-1), anglegrav(2:imax-1) ./ pi .* 180, '-b');
    set(Langle(1), 'LineWidth', 2)
    set(Langle(2), 'LineWidth', 1)
    set(Langle(3), 'LineWidth',1)
end
xlim([min(iters), iters(imax)]);
ylim([0 90]);
grid on
% set(Lmdifn, 'LineWidth', 1)
title('angle between consecutive kernels (degrees)');

% step size
subplot(9,2,[15 17])
istepmax = min(length(step),imax-1);
Lstep = semilogy(iters(1:istepmax), step(1:istepmax),'k');
xlim([min(iters), iters(imax)]);
grid on
set(Lstep, 'LineWidth', 1)
title('step size');

% label for all plots, really, though technically belonging to Lstep
xlabel('iteration no.');

%% - right column

% misfit development
subplot(9,2,[2,6])
Rmisfit = semilogy(stepcumsum(1:imax), misfit(1:imax), 'k',... 
        stepcumsum(1:imax), misfitseis(1:imax), '-r', ...
        stepcumsum(1:imax), misfitgrav(1:imax), '-b');
    set(Rmisfit(1), 'LineWidth', 2)
    set(Rmisfit(2), 'LineWidth', 1)
    set(Rmisfit(3), 'LineWidth',1)
%     h(1).LineWidth = 2;
xlim([min(stepcumsum), stepcumsum(imax)]);
grid on
xlabel('distance travelled in misfit landscape');
title('misfit (normalised)');
% legend('total (used) misfit', 'seismic misfit', 'gravity misfit', 'Location', 'Northeast');

% angle between consecutive descent directions vs. path
subplot(9,2,[12 14])
% IKS = [stepcumsum(2:imax-1);stepcumsum(2:imax-1);stepcumsum(2:imax-1)]';
% EI = [InvProps.Ktotaal_angle(2:imax-1);InvProps.Kseis_angle(2:imax-1); ...
%       InvProps.Kg_angle(2:imax-1)]';
%   Rangle = stem(IKS, EI);
if imax > 2
Rangle = plot(stepcumsum(2:imax-1), angletot(2:imax-1) ./ pi .* 180, '--ko',...
              stepcumsum(2:imax-1), angleseis(2:imax-1) ./ pi .* 180, '--ro', ...
              stepcumsum(2:imax-1), anglegrav(2:imax-1) ./ pi .* 180, '--bo');
    set(Rangle(1), 'LineWidth', 2)
    set(Rangle(2), 'LineWidth', 1)
    set(Rangle(3), 'LineWidth',1)
end
xlim([min(stepcumsum), stepcumsum(imax)]);
ylim([0 90]);
grid on
% set(Lmdifn, 'LineWidth', 1)
title('angle between consecutive kernels (degrees)');
xlabel('distance travelled in misfit landscape');
% yt=get(gca,'ytick')
% for k=1:numel(yt);
% yt1{k}=sprintf('%4.1f°',yt(k));
% end
% set(gca,'yticklabel',yt1);

% % step size
% subplot(4,2,8)
% Lstep = plot(step,'k');
% xlim([min(iters), iters(imax)]);
% set(Lstep, 'LineWidth', 1)
% title('step size');

% label for all plots, really, though technically belonging to Lstep



end