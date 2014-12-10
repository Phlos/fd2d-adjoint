function fig_inv = plot_inversion_development(misfit, misfit_seis, misfit_g, step, modeldifn, Model_end, Model_start, Model_real , middle)

%% prep

set_figure_properties_bothmachines;

% input_parameters;
iters = 1:(length(misfit));
niter = length(iters);
for i = 1:niter
misfitseis1(i) = misfit_seis(i).normd;
end
for i = 1:niter
misfitgrav1(i) = misfit_g(i).normd;
end


misfitseis = misfitseis1 ./ misfitseis1(1);
misfitgrav = misfitgrav1 ./ misfitgrav1(1);
if ( ~(misfit(1) == 1) && ~isequal(misfit, misfitseis1) && ~isequal(misfit, misfitgrav1) )
misfit = misfitseis + misfitgrav;
else
    misfit = misfit ./ misfit(1);
end

% if misfit(1) == 2
%     disp 'the used misfit was normalised to gravity and seis'
%     misfit = misfit ./ (0.5* misfit(1));
%     misfitseis = misfitseis ./ misfitseis(1);
%     misfitgrav = misfitgrav ./ misfitgrav(1);
% else
%     disp 'the used misfit was neither seis nor grav'
%     misfit = misfit ./ misfit(1);
%     misfitseis = misfitseis ./ misfitseis(1);
%     misfitgrav = misfitgrav ./ misfitgrav(1);
% end

% if(misfitseis == misfit)
%     disp 'the used misfit was the seismic misfit'
%     misfit = misfit ./ misfit(1);
%     misfitseis = misfitseis ./ misfitseis(1);
%     misfitgrav = misfitgrav ./ misfitgrav(1);
% elseif (misfitgrav == misfit)
%     disp 'the used misfit was the gravity misfit... WEIRD?!!'
%     misfit = misfit ./ misfit(1);
%     misfitseis = misfitseis ./ misfitseis(1);
%     misfitgrav = misfitgrav ./ misfitgrav(1);
% elseif misfit(1) == 2
%     disp 'the used misfit was normalised to gravity and seis'
%     misfit = misfit ./ (0.5* misfit(1));
%     misfitseis = misfitseis ./ misfitseis(1);
%     misfitgrav = misfitgrav ./ misfitgrav(1);
% elseif misfit(1) == 2
% else
%     disp 'the used misfit was neither seis nor grav'
% end

%% make alllll the nice plots in one effort

fig_inv = figure;
pos_invplot =  [ 1082         474        500        1046];
set(fig_inv, 'OuterPosition', pos_invplot)

%% - left column: misfit development, step length dev., model norm dev.

% misfit development
% subplot(4,4,[1,5])
subplot(4,1,[1,2])
% semilogy(iters, misfit ./ misfit(1), 'k',... 'Linewidth', 2, ...
h = semilogy(iters, misfit, 'k',... 'Linewidth', 2, ...
        iters, misfitseis, '-r', ...
        iters, misfitgrav, '-b');
    set(h(1), 'LineWidth', 2)
    set(h(2), 'LineWidth', 1)
    set(h(3), 'LineWidth',1)
%     h(1).LineWidth = 2;
% xlabel('iteration no.');
% ylabel('misfit (normalised)');
title('misfit (normalised)');
legend('total (used) misfit', 'seismic misfit', 'gravity misfit', 'Location', 'Northoutside');

% step size
% subplot(4,4,9)
subplot(4,1,3)
h = plot(step,'k');
set(h, 'LineWidth', 1)

% xlabel('iteration no.');
title('step size');

% norm of |model - prev model|
% subplot(4,4,13)
subplot(4,1,4)
% plot(modeldifn,'k');
h = semilogy(modeldifn,'k');
set(h, 'LineWidth', 1)
title('|current - previous model| (normalised)');
xlabel('iteration no.');

% %% right column: the models
% 
% % preparing subplots
% % for i = 1:3
% %     a(i) = subplot(3,4,(i-1)+[2,4]);
% % end
% ax(1) = subplot(3,4,[2,4]);
% ax(2) = subplot(3,4,4+[2,4]);
% ax(3) = subplot(3,4,8+[2,4]);
% 
% % getting the real model & plotting it in the right position
% plot_model(Model_real,middle,'rhovsvp');
% h = get(gcf,'Children');
% % put real model into totalfigure
% newh = copyobj(h,fig_inv);
% i=1;
% for j = 1:length(newh)
%     posnewh = get(newh(j),'Position')
%     possub  = get(ax(i),'Position')
%     set(newh(j),'Position',...
%         [posnewh(1) possub(2) posnewh(3) possub(4)])
%     
% end
% 
% delete(ax(i));

end