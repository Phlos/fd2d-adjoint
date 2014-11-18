function [fig_kcomp, Kyxdiff] = plot_different_kernel_filters(X, Z, Knof)

% plot the kernels with differnet filter types

fig_kcomp = figure;
set(fig_kcomp, 'OuterPosition',[1153, 568, 1837, 875]);

% 1
% original kernel
subplot(2,3,1);
pcolor(X,Z,Knof')
shading interp
axis image
colorbar
scale = prctile(abs(Knof(:)), 98);
caxis([-scale scale]);
title('original kernel');

% 2
% filtered Y first
Kfy = filter_kernels(X,Z,Knof',15*334);
subplot(2,3,2);
pcolor(X,Z,Kfy)
shading interp
axis image
colorbar
caxis([-scale scale]);
title('kernel filtered Y first');

% 3
% filtered X first
Kfx = filter_kernels(Z',X',Knof,15*334);
subplot(2,3,3);
pcolor(X,Z,Kfx')
shading interp
axis image
colorbar
caxis([-scale scale]);
title('kernel filtered X first');

% 4
% diff plot between Y first and X first
Kyxdiff = Kfy' - Kfx;
subplot(2,3,4);
% Kginv_filt = filter_kernels(Z',X',Kg,15*334);
pcolor(X,Z,Kyxdiff')
shading interp
axis image
colorbar
caxis([-scale scale]);
title('difference between Y first and X first');

% 5
% filtered with the original filter
filt = fspecial('gaussian',[15 15], 9);
Kforig = conv2(Knof, filt, 'same');
subplot(2,3,5);
pcolor(X,Z,Kforig')
shading interp
axis image
colorbar
% scale = prctile(abs(Kforig(:)), 98);
caxis([-scale scale]);
title('kernel filtered w previous filter');

%6
% diff plot between Xfirst and orig filter
Kyorigdiff = Kfx - Kforig;
subplot(2,3,6);
pcolor(X,Z,Kyorigdiff')
shading interp
axis image
colorbar
caxis([-scale scale]);
title('difference between X first and previous filter');
