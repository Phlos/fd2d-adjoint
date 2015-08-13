function fig_invres = plot_inversion_result_InvTlbx(InvProps, imax)

% plotting the results of the inversion in a good format.
% 3 x 2 subplots with:
% (left column:)
% - misfit(s) vs. iter
% - modeldifnorm vs. iter
% (right column:)
% - misfit vs. cumulative steplength
% - angle vs. misfit path? or iter?!

%% prepare

input_parameters;
set_figure_properties_bothmachines;

% variables

% x axis:
% step = InvProps.step;

% misfits
misfit = InvProps.misfit;
misfitseis = InvProps.misfitseis;
misfitgrav = InvProps.misfitgrav;

% norm of model differences
modeldifn = InvProps.modeldifn;
if(isfield(InvProps,'modeldifnFromTrue'))
    normL2 = InvProps.modeldifnFromTrue;
end

% kernel angles
angletot = InvProps.angle.Ktotal;
if strcmp(use_seis,'yesseis')
    angleseis = InvProps.angle.Kseis;
else
    angleseis = NaN*zeros(size(angletot));
end
if strcmp(use_grav,'yes')
    anglegrav = InvProps.angle.Kg;
else
    anglegrav = NaN*zeros(size(angleseis));
end

% kernel norms
normKseis = InvProps.norm.Kseis;
normKgrav = InvProps.norm.Kg;
normKtot  = InvProps.norm.Ktotal;

% some adaptations
niter = length(misfit); % the actual total number of iters run
imax = min(imax,niter);
iters = 1:niter;

% for ii = 1:imax
%     misfitseis(ii) = misfit_seis{ii}.normd;
%     misfitgrav(ii) = misfit_g(ii).normd;
% end

% stepcumsum = [0, cumsum(abs(step(1:imax-1)))];

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
xlim([min(iters), iters(imax)]);
grid on
legend('all data combined', 'seismic data only', 'gravity data only', ...
       'Location', 'Northeast');
text(0.5, 0.9, 'misfit (normalised)', ...
    'Units', 'normalized', 'HorizontalAlignment','center');

% norm of |model - prev model|
subplot(9,2,[7 9])
if imax > 1
Lmdifn = semilogy(iters(2:imax), modeldifn(2:imax),'k');
end
xlim([min(iters), iters(imax)]);
grid on
set(Lmdifn, 'LineWidth', 1)
text(0.5, 0.9, '|current - previous model| / |1st model|', ...
    'Units', 'normalized', 'HorizontalAlignment','center');


% angle between consecutive kernels - vs. iter
subplot(9,2,[11 13])
if imax > 2
Langle = plot(iters(2:imax), angletot(2:imax) ./ pi .* 180, 'k',... 
        iters(2:imax), angleseis(2:imax) ./ pi .* 180, '-r', ...
        iters(2:imax), anglegrav(2:imax) ./ pi .* 180, '-b');
    set(Langle(1), 'LineWidth', 2)
    set(Langle(2), 'LineWidth', 1)
    set(Langle(3), 'LineWidth',1)
end
xlim([min(iters), iters(imax)]);
ylim([0 90]);
grid on
    text(0.5, 0.9, 'angle between consecutive kernels(\circ)', ...
    'Units', 'normalized', 'HorizontalAlignment','center');

% % step size
% subplot(9,2,[15 17])
% istepmax = min(length(step),imax-1);
% Lstep = semilogy(iters(1:istepmax), step(1:istepmax),'k');
% xlim([min(iters), iters(imax)]);
% grid on
% set(Lstep, 'LineWidth', 1)
%     text(0.5, 0.9, 'step size', ...
%     'Units', 'normalized', 'HorizontalAlignment','center');

% label for all plots, really, though technically belonging to Lstep
xlabel('iteration no.');

%% - right column: vs cumulative steplength

% % misfit development
% subplot(9,2,[2,6])
% Rmisfit = semilogy(stepcumsum(1:imax), misfit(1:imax), 'k',... 
%         stepcumsum(1:imax), misfitseis(1:imax), '-r', ...
%         stepcumsum(1:imax), misfitgrav(1:imax), '-b');
%     set(Rmisfit(1), 'LineWidth', 2)
%     set(Rmisfit(2), 'LineWidth', 1)
%     set(Rmisfit(3), 'LineWidth',1)
% %     h(1).LineWidth = 2;
% xlim([min(stepcumsum), stepcumsum(imax)]);
% grid on
% % xlabel('distance travelled in misfit landscape');
% % title('misfit (normalised)');
%     text(0.5, 0.9, 'misfit (normalised)', ...
%     'Units', 'normalized', 'HorizontalAlignment','center');
% % legend('total (used) misfit', 'seismic misfit', 'gravity misfit', 'Location', 'Northeast');

% model diff with real model vs. step length
if(exist('normL2','var'))
    subplot(9,2,[8 10]);
    LmdifReal = semilogy(iters(1:imax), normL2(1:imax),'k');
    xlim([min(iters) iters(imax)]);
    grid on
    set(LmdifReal, 'LineWidth', 1)
    text(0.5, 0.9, '|current - real| / |real - bg|', ...
    'Units', 'normalized', 'HorizontalAlignment','center');
%     title('|current - real| / |real - bg|');
end

% % angle between consecutive descent directions vs. cumulative step
% subplot(9,2,[12 14])
% if imax > 2
% Rangle = plot(stepcumsum(2:imax-1), angletot(2:imax-1) ./ pi .* 180, '--ko',...
%               stepcumsum(2:imax-1), angleseis(2:imax-1) ./ pi .* 180, '--ro', ...
%               stepcumsum(2:imax-1), anglegrav(2:imax-1) ./ pi .* 180, '--bo');
%     set(Rangle(1), 'LineWidth', 2)
%     set(Rangle(2), 'LineWidth', 1)
%     set(Rangle(3), 'LineWidth',1)
% end
% xlim([min(stepcumsum), stepcumsum(imax)]);
% ylim([0 180]);
% grid on
% % set(Lmdifn, 'LineWidth', 1)
% text(0.5, 0.9, 'angle between consecutive kernels (\circ)', ...
%     'Units', 'normalized', 'HorizontalAlignment','center');



% kernel magnitude
subplot(9,2,[16 18])
RKnorm = semilogy(iters(1:imax), normKtot(1:imax), '-k',...
              iters(1:imax), normKseis(1:imax), '-r',...
              iters(1:imax), normKgrav(1:imax), '-b');
set(RKnorm(1), 'LineWidth', 2)
set(RKnorm(2), 'LineWidth', 1)
set(RKnorm(3), 'LineWidth',1)
xlim([min(iters), iters(imax)]);
grid on;
text(0.5, 0.9, 'norm of kernels', ...
    'Units', 'normalized', 'HorizontalAlignment','center');
% title('norm of kernels');

xlabel('distance travelled in misfit landscape');




end