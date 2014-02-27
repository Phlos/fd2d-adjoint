% function plot_model(fig_handle,X,Z,mu,lambda,rho)

load 'propagation/cm_model.mat';

fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod) 
set(gca,'FontSize',14)

subplot(1,3,1)
pcolor(X,Z,mu');
colormap(cm_model);
axis image
shading flat
title('mu [N/m^2]','FontSize',14)
xlabel('x [m]','FontSize',14);
ylabel('z [m]','FontSize',14);
colorbar
    
subplot(1,3,3)
pcolor(X,Z,rho');
colormap(cm_model);
axis image
shading flat
title('rho [kg/m^3]','FontSize',14)
xlabel('x [m]','FontSize',14);
ylabel('z [m]','FontSize',14);
colorbar

subplot(1,3,2)
pcolor(X,Z,lambda');
colormap(cm_model);
axis image
shading flat
title('lambda [N/m^2]','FontSize',14)
xlabel('x [m]','FontSize',14);
ylabel('z [m]','FontSize',14);
colorbar