function plot_kernel(X,Z,kernel,kname,cmaxtype,cmax,stf_PSV)

% plot some random sensitivity kernel
%
% input:
% ------
% X,Z:      X and Z grids
% kernel:   the sensitivity kernel you want to plot
% kname:    kernel name (should be string)
% cmaxtype: 'perc' or 'fixed'. Should the maxima of the kernel plot be
%           determined by a percentile of the kernel values ('perc'), or do
%           you want a fixed value ('fixed')
% cmax:     either the percentile of the colours for the colorbar maximum
%           in case cmaxtype = 'perc', or the maximum value of the
%           colourscale in case cmaxtype = 'fixed'.
% 
% output:
% -------
% a plot of the sensitivity kernel
%==========================================================================

load 'cm_velocity.mat';

if strcmp(cmaxtype,'perc')
    scale = prctile(kernel(:),cmax);
elseif strcmp(cmaxtype,'fixed')
    scale = cmax;
else
    error('Invalid cmaxtype. Allowed values are ''perc'' or ''fixed''.')
end

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


caxis([-scale scale]);

text(0.95*max(X(:)),0.92*max(Z(:)),['max = \pm', num2str(scale,'%3.1e')], ... 
                      'HorizontalAlignment','right')

hold off;

end