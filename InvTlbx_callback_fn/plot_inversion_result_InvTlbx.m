function fig_invres = plot_inversion_result_InvTlbx(InvProps, varargin)

% plotting the results of the inversion in a good format.
% 3 x 2 subplots with:
% (left column:)
% - misfit(s) vs. iter
% - modeldifnorm vs. iter
% (right column:)
% - misfit vs. cumulative steplength
% - angle vs. misfit path? or iter?!

%% prepare

if numel(varargin) == 0
    imax = numel(InvProps.misfit);
elseif numel(varargin) == 1
    imax = varargin{1};
else
    error('Panic!');
end
    

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
if(isfield(InvProps,'L2norm_normd'))
    normL2.rvv = InvProps.L2norm_normd.total_rvv;
    if isfield(InvProps.L2norm_normd, 'total_rml')
        normL2.rml = InvProps.L2norm_normd.total_rml;
    else
        normL2.rml = (InvProps.L2norm_normd.rho + InvProps.L2norm_normd.mu + InvProps.L2norm_normd.lambda) / 3;
    end
    normL2Par.rho    = InvProps.L2norm_normd.rho;
    normL2Par.mu     = InvProps.L2norm_normd.mu;
    normL2Par.lambda = InvProps.L2norm_normd.lambda;
    normL2Par.vs     = InvProps.L2norm_normd.vs;
    normL2Par.vp     = InvProps.L2norm_normd.vp;
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
if (feature('showfigurewindows') == 0)
    set(fig_invres, 'PaperUnits', 'points');
    set(fig_invres, 'PaperPosition', pos_invres);
end

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
if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
grid on
legenda = legend('all data combined', 'seismic data only', 'gravity data only');
lp = get(legenda, 'Position');
set(legenda, 'Position', [lp(1)+0.2 lp(2) lp(3), lp(4)]);

text(0.5, 0.9, 'misfit (normalised)', ...
    'Units', 'normalized', 'HorizontalAlignment','center');

% norm of |model - prev model|
subplot(9,2,[7 9])
if imax > 1
Lmdifn = semilogy(iters(2:imax), modeldifn(2:imax),'k');
end
xlim([min(iters), iters(imax)]);
if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
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
if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
ax = gca;
ax.YTick = 0:30:180;
ylim([0 180]);
grid on
    text(0.5, 0.9, 'angle between consecutive kernels(\circ)', ...
    'Units', 'normalized', 'HorizontalAlignment','center');


% kernel magnitude
subplot(9,2,[15 17])
RKnorm = semilogy(iters(1:imax), normKtot(1:imax), '-k',...
              iters(1:imax), normKseis(1:imax), '-r',...
              iters(1:imax), normKgrav(1:imax), '-b');
set(RKnorm(1), 'LineWidth', 2)
set(RKnorm(2), 'LineWidth', 1)
set(RKnorm(3), 'LineWidth',1)
xlim([min(iters), iters(imax)]);
if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
grid on;
text(0.5, 0.9, 'norm of kernels', ...
    'Units', 'normalized', 'HorizontalAlignment','center');
% title('norm of kernels');

% label for all plots, really, though technically belonging to Lstep
xlabel('iteration no.');

%% - right column: other stuff

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

% total model diff with real model development vs iteration
if(exist('normL2','var'))
    subplot(9,2,[8 10]); box on;
    hold on;
    RmdifReal_rvv = plot(iters(1:imax), normL2.rvv(1:imax));
    RmdifReal_rml = plot(iters(1:imax), normL2.rml(1:imax));
    eilim = ylim;
    if eilim(2) > 10
	    eilim(2) = 2;
    end
    ylim([0 max(1, eilim(2))]);
    xlim([min(iters) iters(imax)]);
    if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
    grid on
    set(RmdifReal_rvv, 'LineWidth', 1)
    set(RmdifReal_rml, 'LineWidth', 1)
    text(0.5, 0.9, '| current - real | / | real - bg |', ...
    'Units', 'normalized', 'HorizontalAlignment','center');
    legenda = legend('mean  rvv ', 'mean  rml ');
    
    lp = get(legenda, 'Position');
    set(legenda, 'Position', [lp(1) lp(2)+0.075 lp(3), lp(4)]);
%     title('|current - real| / |real - bg|');
end

% % L2 norm between current model and true model per parameter
if(exist('normL2Par', 'var'))
    subplot(9,2,[12 14]); box on;
    hold on;
    RdifrealPar_rho = plot(iters(1:imax), normL2Par.rho(1:imax));
%     if sum(normL2Par.vs) > 1e-5*sum(normL2Par.rho)
        RdifrealPar_vs = plot(iters(1:imax), normL2Par.vs(1:imax));
        RdifrealPar_vp = plot(iters(1:imax), normL2Par.vp(1:imax));
        RdifrealPar_mu = plot(iters(1:imax), normL2Par.mu(1:imax));
        RdifrealPar_lm = plot(iters(1:imax), normL2Par.lambda(1:imax));
%     end
    eilim = ylim;
    if eilim(2) > 10
	    eilim(2) = 5;
    end
    ylim([0 max(1, eilim(2))]);
    xlim([min(iters) iters(imax)]);
    if imax>80; set(gca, 'XTick', 0:20:iters(imax)); end;
    eilim = ylim;
    if eilim(1) < 5e-2
        ylim([5e-2 eilim(2)]);
    end
    grid on;
    text(0.5, 0.9, '| current - real | / | real - bg |', ...
       'Units', 'normalized', 'HorizontalAlignment','center');
   
    legenda = legend([RdifrealPar_rho, RdifrealPar_vs, RdifrealPar_vp, RdifrealPar_mu, RdifrealPar_lm], ...
        'density', 'S velocity', 'P velocity', 'mu', 'lambda'); %, ...
%         'Location', 'southoutside', 'orientation', 'horizontal');
    lp = get(legenda, 'Position');
    set(legenda, 'Position', [lp(1) lp(2)-0.2 lp(3), lp(4)]);
end


xlabel('iteration no.');


% % step size
% subplot(9,2,[16 18])
% istepmax = min(length(step),imax-1);
% Lstep = semilogy(iters(1:istepmax), step(1:istepmax),'k');
% xlim([min(iters), iters(imax)]);
% grid on
% set(Lstep, 'LineWidth', 1)
%     text(0.5, 0.9, 'step size', ...
%     'Units', 'normalized', 'HorizontalAlignment','center');






end
