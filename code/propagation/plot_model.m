% function plot_model(fig_handle,X,Z,mu,lambda,rho)

load 'propagation/cm_model.mat';

fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod) 
% set(gca,'FontSize',14)

subplot(1,3,1)

hold on
pcolor(X,Z,mu');

%- colour scale
if all(mu == mu(1))
    cmax = mu(1) + 0.01*mu(1);
    cmin = mu(1) - 0.01*mu(1);
else
    cmax = max(mu(:));
    cmid = mode(mu(:));
    cmin = 2*cmid - cmax;
end
caxis([cmin cmax]);


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


subplot(1,3,2)
hold on
pcolor(X,Z,lambda');

%- colour scale
if all(lambda == lambda(1))
    cmax = lambda(1) + 0.01*lambda(1);
    cmin = lambda(1) - 0.01*lambda(1);
else
    cmax = max(lambda(:));
    cmid = mode(lambda(:));
    cmin = 2*cmid - cmax;
end
caxis([cmin cmax]);

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




subplot(1,3,3)
hold on
pcolor(X,Z,rho');

%- colour scale
if all(rho == rho(1))
    cmax = rho(1) + 0.01*rho(1);
    cmin = rho(1) - 0.01*rho(1);
else
    cmax = max(rho(:));
    cmid = mode(rho(:));
    cmin = 2*cmid - cmax;
end
caxis([cmin cmax]);

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

