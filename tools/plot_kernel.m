function plot_kernel(X,Z,kernel,kname,prc)

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


% cmax = min(max(kernel(:)),8*std(kernel(:)));
% caxis([-cmax cmax]);

% experimenting with the caxis (3) -- location of maximum value
% could still use it like in experiment (1) blanking out the max region
% [value, k] = max(kernel(:));
% [ie, jee] = ind2sub(size(kernel), k);
% max_loc = [X(1,ie),Z(jee)];

% K_no_origsrc = kernel;
% blank=min(size(X))/15;
% K_no_origsrc( ie-blank:ie+blank , jee-blank:jee+blank ) = 0; 

% ALTERNATIVE:
% the same but for every source instead of only the maximum location
% for i=1:ns_orig
%     iks = orig_x_id(i);
%     zet = orig_z_id(i);
%     K_no_origscr( iks-blank:iks+blank , zet-blank:zet+blank ) = 0;
% end


% figure(figname);
pcolor(X,Z,kernel');
% if(test);
%     pcolor(X,Z,K_no_origsrc');
% end
title({['kernel for ', kname]; '(source direction x,z = (1,-5) )'});
shading interp
colorbar
axis image
% colormap(flipud(cm));
colormap(cm);
hold on;


% p = plot(X(1,ie),Z(jee),'kx');
% set(p,'Color','green');

% cmax = max(max(abs(K_no_origsrc)));
cmax = prctile(kernel(:),prc);
caxis([-cmax cmax]);

hold off;

end