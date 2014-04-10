% function plot_model(fig_handle,X,Z,mu,lambda,rho)

load 'propagation/cm_model.mat';

fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod) 
% set(gca,'FontSize',14)

subplot(1,3,1)

hold on
pcolor(X,Z,mu');
        for k=1:length(src_x)
            plot(src_x(k),src_z(k),'kx')
        end
        
        for k=1:length(rec_x)
            plot(rec_x(k),rec_z(k),'ko')
        end
colormap(cm_model);
axis image
shading flat
title('mu [N/m^2]')
xlabel('x [m]');
ylabel('z [m]');
colorbar
hold off;

subplot(1,3,3)
hold on
pcolor(X,Z,rho');
        for k=1:length(src_x)
            plot(src_x(k),src_z(k),'kx')
        end
        
        for k=1:length(rec_x)
            plot(rec_x(k),rec_z(k),'ko')
        end
colormap(cm_model);
axis image
shading flat
title('rho [kg/m^3]')
xlabel('x [m]');
ylabel('z [m]');
colorbar
hold off;

subplot(1,3,2)
hold on
pcolor(X,Z,lambda');
        for k=1:length(src_x)
            plot(src_x(k),src_z(k),'kx')
        end
        
        for k=1:length(rec_x)
            plot(rec_x(k),rec_z(k),'ko')
        end
colormap(cm_model);
axis image
shading flat
title('lambda [N/m^2]')
xlabel('x [m]');
ylabel('z [m]');
colorbar
hold off;

