function plot_kernel(X,Z,kernel,kname,cmaxtype,cmax,varargin)

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

input_parameters;
checkargs(varargin);


if strcmp(cmaxtype,'perc')
    scale = prctile(abs(kernel(:)),cmax);
elseif strcmp(cmaxtype,'fixed')
    scale = cmax;
else
    error('Invalid cmaxtype. Allowed values are ''perc'' or ''fixed''.')
end

% figure(figname);
pcolor(X,Z,kernel');

% title({['kernel for ', kname]; ...
%         ['(source direction x,z = (',num2str(stf_PSV), ') )']});
title({['kernel for ', kname]});
shading interp
% colorbar
axis image
% colormap(flipud(cm));
% colormap(cm);
hold on;

for k=1:length(src_x)
    plot(src_x(k),src_z(k),'kx','LineWidth',0.3,'MarkerSize',4)
end

for k=1:length(rec_x)
    plot(rec_x(k),rec_z(k),'ko','LineWidth',0.3,'MarkerSize',4)
end

caxis([-scale scale]);

% text(0.95*max(X(:)),0.92*max(Z(:)),['max = \pm', num2str(scale,'%3.1e')], ... 
%                       'HorizontalAlignment','right', ...
%                       'Color',[0.5 0.5 0.5])
text(0.95,0.92,['max = \pm', num2str(scale,'%3.1e')], ... 
                      'Units', 'normalized', ...
                      'HorizontalAlignment','right', ...
                      'Color',[0.5 0.5 0.5])

hold off;

end

function checkargs(arg)

narg = length(arg);
if narg == 1
    colormap(arg{1})
elseif narg == 0;
    load './code/propagation/cm_velocity.mat';
    colormap(cm);
else
    error('wrong number of variable arguments')
end

end