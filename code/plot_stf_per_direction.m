% plot all the source time functions

figure;
for i = 1:3
    subplot(3,1,i)
    hold on
    for j= 1:ns
        plot(t,squeeze(stf(i,j,:)));
    end
end