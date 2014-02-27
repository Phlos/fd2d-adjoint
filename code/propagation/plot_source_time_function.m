% plots the source time function in the format described here.

fig_stf = figure;
set(fig_stf,'OuterPosition',pos_stf);
set(gca,'FontSize',20)
hold on

plot(t,stf,'k');
    xlabel('time [s]','FontSize',20);
    title('source-time function','FontSize',20);