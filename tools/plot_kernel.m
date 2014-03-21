function plot_kernel(X,Z,kernel,kname,prc,stf_PSV)

% plot some random sensitivity kernel
%
% input:
% ------
% X,Z:      X and Z grids
% kernel:   the sensitivity kernel you want to plot
% kname:    kernel name (should be string)
% prc:      percentile of the colours for the colorbar maximum
% 
% output:
% -------
% a plot of the sensitivity kernel
%==========================================================================

load 'cm_velocity.mat';



% figure(figname);
pcolor(X,Z,kernel');

title({['kernel for ', kname]; ...
        ['(source direction x,z = (',num2str(stf_PSV), ') )']});
shading interp
% colorbar
axis image
% colormap(flipud(cm));
colormap(cm);
hold on;


% p = plot(X(1,ie),Z(jee),'kx');
% set(p,'Color','green');

% cmax = max(max(abs(K_no_origsrc)));
cmax = prctile(kernel(:),prc);

caxis([-cmax cmax]);

text(0.95*max(X(:)),0.92*max(Z(:)),['max = \pm', num2str(cmax,'%3.1e')], ... 
                      'HorizontalAlignment','right')

hold off;

end