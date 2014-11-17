% plot the original model, the current test model, and the kernels

% Input:
%
% - original model params
% - current test model params
% - sensitivity kernels

% OUTPUT:
%
% - a plot with real model, test model and sensitivity kernels

%==========================================================================

function plot_realmodel_testmodel_kernel(realmodelnr,testmodelnr,K)

% initialise some stuff
path(path,'../code')
load 'propagation/cm_model.mat';
cm = cm_model;

input_parameters;
[X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);
set_figure_properties_bothmachines;


fig_rtk = figure;
set(fig_rtk,'OuterPosition',pos_rtk);

%- plot original model (top row)

[mu,rho,lambda]=define_material_parameters(nx,nz,realmodelnr);

for i = 1:3
    
    if i==1
        parameter = mu;
        titulo = '\mu [N/m^2]';
    elseif i==2
        parameter = lambda;
        titulo = '\lambda [N/m^2]';
    elseif i==3
        parameter = rho;
        titulo = '\rho [kg/m^3]';
    end
    
    subplot(3,3,i)
    
    hold on
    pcolor(X,Z,parameter');
    
    %- colour scale
    if all(parameter == parameter(1))
        cmax = parameter(1) + 0.01*parameter(1);
        cmin = parameter(1) - 0.01*parameter(1);
    else
        cmax = max(parameter(:));
        cmid = mode(parameter(:));
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
    
    %- titles & labels
    title(titulo);
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    hold off;
    
    
    
end




%- plot the test model (middle row)
[mu,rho,lambda]=define_material_parameters(nx,nz,testmodelnr);

for i = 1:3
    
    if i==1
        parameter = mu;
        titulo = '\mu [N/m^2]';
    elseif i==2
        parameter = lambda;
        titulo = '\lambda [N/m^2]';
    elseif i==3
        parameter = rho;
        titulo = '\rho [kg/m^3]';
    end
    
    subplot(3,3,3+i)
    
    hold on
    pcolor(X,Z,parameter');
    
    %- colour scale
    if all(parameter == parameter(1))
        cmax = parameter(1) + 0.01*parameter(1);
        cmin = parameter(1) - 0.01*parameter(1);
    else
        cmax = max(parameter(:));
        cmid = mode(parameter(:));
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
    
    %- titles & labels
    title(titulo);
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    hold off;
    
    
end

%- plot the (total) sensitivity kernels (bottom row)

% for parameter = {'rho' 'mu' 'lambda'};

for i = 1:3
    
    if i==1
        parameter = {'mu'};
        titulo = '\mu [N/m^2]';
    elseif i==2
        parameter = {'lambda'};
        titulo = '\lambda [N/m^2]';
    elseif i==3
        parameter = {'rho'};
        titulo = '\rho [kg/m^3]';
    end
    
    subplot(3,3,6+i)
    plot_kernel(X,Z,K.(parameter{1}).total,titulo,'perc',99.97,stf_PSV,cm_model);
    
    colorbar
    
end

end