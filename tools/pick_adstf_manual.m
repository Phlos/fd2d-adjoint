function [left, right] = pick_adstf_manual(t, v_rec, v_obs, n, fig_manual)

fprintf(1,'station number %d\n',n)

figure(fig_manual);
clf;
subplot(3,1,1)
plot(t,v_rec(n,:),'k')
hold on
plot(t,v_obs(n,:),'r')
plot(t,v_rec(n,:)-v_obs(n,:),'k--','LineWidth',2)
% hold off

title(['receiver ' num2str(n) ' ,synth - black, obs - red, diff - dashed'])
xlabel('t [s]')
ylabel('velocity or displacement')

disp('select left window');
[left,~]=ginput(1);
disp('select right window');
[right,~]=ginput(1);

% hold off;
ei_as = get(gca,'ylim');
plot([left, left], ei_as, '--');
plot([right, right], ei_as, '--');
hold off

end