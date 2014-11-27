function fig_knl = plot_kernels_rho_vs_vp(K)

% plot kernels for SH and P-SV in rho, mu, lambda

%- initialise input -------------------------------------------------------
path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
[mu,rho,lambda]=define_material_parameters(nx,nz,11);

%- plot -------------------------------------------------------------------

% figure position and size
fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)

% P-SV
subplot(3,2,1)
plot_kernel(X,Z,K.rho2.PSV,'rho - P-SV','perc',99.95,stf_PSV);
subplot(3,2,3)
plot_kernel(X,Z,K.vs2.PSV,'vs - P-SV','perc',99.95,stf_PSV);
subplot(3,2,5)
plot_kernel(X,Z,K.vp2.PSV,'vp - P-SV','perc',99.95,stf_PSV);
% SH
subplot(3,2,2)
plot_kernel(X,Z,K.rho2.SH,'rho - SH','perc',99.97,stf_PSV);
subplot(3,2,4)
plot_kernel(X,Z,K.vs2.SH,'v_s - SH','perc',99.98,stf_PSV);
subplot(3,2,6)
plot_kernel(X,Z,zeros(size(X')),'v_p - SH','perc',100,stf_PSV);

% colourbar
h = colorbar('horiz');
set(h,'XTick',[-1;0;1])
set(h,'Xticklabel',{'-max','0', 'max'})
set(h,'Position',[0.3 0.05 0.4 0.02])

end