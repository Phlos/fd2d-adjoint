% plot the test values of stf, DSX, vx, DSX before stf is added
%
% ultimate goal of this plot is to find out why vx is so much smaller in my
% model runs than in those of Tromp et al.

figure;
plot(t,test.stf,t,test.DSX,t,test.vx,t,test.DSX_before_stf);
max(test.vx)