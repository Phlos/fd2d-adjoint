% plot kernels for SH and P-SV in rho, mu, lambda

function plot_kernels_rho_mu_lambda(K)

path(path,'../code');
path(path,'../input');

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
% [mu,rho,lambda]=define_material_parameters(nx,nz,11);
set_figure_properties_bothmachines;



fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)

% P-SV
subplot(3,2,1)
plot_kernel(X,Z,K.rho.PSV,'rho - P-SV','perc',99.95,stf_PSV);
subplot(3,2,3)
plot_kernel(X,Z,K.mu.PSV,'mu - P-SV','perc',99.95,stf_PSV);
subplot(3,2,5)
plot_kernel(X,Z,K.lambda.PSV,'lambda - P-SV','perc',99.95,stf_PSV);
% SH
subplot(3,2,2)
plot_kernel(X,Z,K.rho.SH,'rho - SH','perc',99.95,stf_PSV);
subplot(3,2,4)
plot_kernel(X,Z,K.mu.SH,'mu - SH','perc',99.95,stf_PSV);
subplot(3,2,6)
plot_kernel(X,Z,zeros(size(X')),'lambda - SH','perc',100,stf_PSV);

h = colorbar('horiz');
set(h,'XTick',[-1;0;1])
set(h,'Xticklabel',{'-max','0', 'max'})
set(h,'Position',[0.3 0.05 0.4 0.02])

end