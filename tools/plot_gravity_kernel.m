function fig_Kg = plot_gravity_kernel(Kg)

% input
input_parameters;
[X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);

load './code/propagation/cm_velocity.mat';

% smooth kernel
Kg_sm = filter_2Dfield(Kg, smoothgwid);

%- output figure
fig_Kg = figure;
pcolor(X,Z, Kg_sm')
colormap(cm);
shading interp
axis image
titel = 'gravity kernel';
title(titel);
colorbar
scale = prctile(abs(Kg_sm(:)),98);
caxis([-scale scale]);

end