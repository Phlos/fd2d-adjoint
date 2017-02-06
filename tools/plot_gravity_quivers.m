% plot gravity field nicely
function fig_grav = plot_gravity_quivers(rec_g, g1, g2, X, Z, rho)
%
% function to plot gravity quivers of the difference g1 - g2
%
% fig_grav = plot_gravity_quivers(rec_g, g1, g2, X, Z, rho)

%% preparation
load './code/propagation/cm_model.mat';

fig_grav = figure;
clf;

%% actual plotting

hold on;

% then plot image
X = X ./ 1000;
Z = Z ./ 1000;
rec_g.x = rec_g.x ./1000;
rec_g.z = rec_g.z ./1000;
pcolor(X,Z,rho');
colormap(cm_model);
shading flat; axis image;
xlabel('x [km]')
ylabel('z [km]')
colorbar;

% plot the gravity receivers
plot(rec_g.x, rec_g.z, 'ko');

% arrow plots ('quiver') with gravity info
% h1 = quiver(rec_g.x,rec_g.z,g1.x, g1.z, 3,'color','k'); % the '0' turns off auto-scaling
% h2 = quiver(rec_g.x,rec_g.z,g2.x, g2.z, 3,'color','r');
difference.x = (g1.x - g2.x);
difference.z = (g1.z - g2.z);
% h3 = quiver(rec_g.x,rec_g.z, difference.x, difference.z, 0,'color','g');
h3 = quivers(rec_g.x,rec_g.z,difference.x, difference.z, 0.5, 3, 'm/s^2','g');

% scaling of quivers
maxi.values = max([g1.x(:);g1.z(:)]);
maxi.domain = max([X(:); Z(:)]);
maxi.toti   = maxi.domain / maxi.values;
qscale = 10*maxi.toti ; % scaling factor for all vectors

% This doesn't work for some reason.
% hendels = {'h1'; 'h2'; 'h3'};
% for i=1:length(hendels)
%     hU = get(hendels{i},'UData') ;
%     hV = get(hendels{i},'VData') ;
%     set(hendels{i},'UData',qscale*hU,'VData',qscale*hV)
% end

% hU = get(h1,'UData') ;
% hV = get(h1,'VData') ;
% set(h1,'UData',qscale*hU,'VData',qscale*hV)
% hU = get(h2,'UData') ;
% hV = get(h2,'VData') ;
% set(h2,'UData',qscale*hU,'VData',qscale*hV)
% hU = get(h3,'UData') ;
% hV = get(h3,'VData') ;
% set(h3,'UData',qscale*hU,'VData',qscale*hV)

end