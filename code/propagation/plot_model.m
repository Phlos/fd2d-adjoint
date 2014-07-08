function fig_mod = plot_model(Model)

% function that plots the model in rho mu lambda parametrisation.
% NOTE: Currenly this is a rather ugly function and I could make it much
% nicer by putting the plotting commands within a for loop over rho mu
% lambda.
%
% INPUT:
% Model:    struct containing .rho .mu .lambda
%
% OUTPUT:   
% - figure with the model plotted. Colour scale are the actual max and min
%   of the parameter values, not some standard deviation.

% format long

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
set_figure_properties_maggi;
mu = Model.mu;
rho = Model.rho;
lambda = Model.lambda;

load 'propagation/cm_model.mat';

fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod) 
% set(gca,'FontSize',14)

subplot(1,3,1)

hold on
pcolor(X,Z,mu');

% difference_mumax_mumin = max(mu(:)) - min(mu(:));

%- colour scale
if all(mu == mu(1))
    cmax = mu(1) + 0.01*mu(1);
    cmin = mu(1) - 0.01*mu(1);
%     disp 'bips!!! all mu are the same!'
else
    % max and min are calculated in this way so that the most common value
    % (i.e. the background value) is white, and that the extreme coulours
    % are determined by whichever of max and min is farthest off from the
    % background colour.
    cmid = mode(mu(:));
    cmax = cmid + max(abs(mu(:) - cmid));
    cmin = cmid - max(abs(mu(:) - cmid));
%     cmax = max(mu(:));
%     cmin = 2*cmid - cmax;
%     disp 'the mu are not all the same'
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
    % max and min are calculated in this way so that the most common value
    % (i.e. the background value) is white, and that the extreme coulours
    % are determined by whichever of max and min is farthest off from the
    % background colour.
    cmid = mode(lambda(:));
    cmax = cmid + max(abs(lambda(:) - cmid));
    cmin = cmid - max(abs(lambda(:) - cmid));
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
    % max and min are calculated in this way so that the most common value
    % (i.e. the background value) is white, and that the extreme coulours
    % are determined by whichever of max and min is farthest off from the
    % background colour.
    cmid = mode(rho(:));
    cmax = cmid + max(abs(rho(:) - cmid));
    cmin = cmid - max(abs(rho(:) - cmid));
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

end
