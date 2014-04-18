function plot_kernel(X,Z,kernel,kname,cmaxtype,cmax,stf_PSV,varargin)

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


caxis([-scale scale]);

text(0.95*max(X(:)),0.92*max(Z(:)),['max = \pm', num2str(scale,'%3.1e')], ... 
                      'HorizontalAlignment','right')

hold off;

end

function checkargs(arg)

narg = length(arg);
if narg == 1
    colormap(arg{1})
elseif narg == 0;
    load '../code/propagation/cm_velocity.mat';
    colormap(cm);
else
    error('wrong number of variable arguments')
end

end