% plot gravity field nicely
function fig_grav = plot_gravity_quivers(rec_g, g, g_homog, X, Z, rho)
%
% function to plot gravity quivers of the difference between g and g_homog
%
% fig_grav = plot_gravity_quivers(rec_g, g, g_homog, X, Z, rho)

fig_grav = figure;
clf;

% first plot the gravity receiver since they're outside of the actual img
plot(rec_g.x, rec_g.z, 'ko');
hold on;

% then plot image
pcolor(X,Z,rho');
shading interp;
axis image;

% arrow plots ('quiver') with gravity info
% h1 = quiver(rec_g.x,rec_g.z,g.x, g.z, 0,'color','g'); % the '0' turns off auto-scaling
% h2 = quiver(rec_g.x,rec_g.z,g_homog.x, g_homog.z, 0,'color','r');
difference.x = (g.x - g_homog.x);
difference.z = (g.z - g_homog.z);
% h3 = quiver(rec_g.x,rec_g.z, difference.x, difference.z, 0,'color','g');
h3 = quivers(rec_g.x,rec_g.z,difference.x, difference.z, 0.5, 3, 'm/s^2','g');

% scaling
maxi.values = max([g.x(:);g.z(:)]);
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