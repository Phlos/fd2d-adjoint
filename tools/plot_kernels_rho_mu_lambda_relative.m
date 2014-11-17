% plot relative kernels rho mu lambda

function fig_knl = plot_kernels_rho_mu_lambda_relative(K_rel)

%- initialise input -------------------------------------------------------
path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
[mu,rho,lambda]=define_material_parameters(nx,nz,11);

%- plot -------------------------------------------------------------------

kernels_together = [K_rel.rho.PSV; K_rel.rho.SH; ...
                    K_rel.mu.PSV; K_rel.mu.SH; ...
                    K_rel.lambda.PSV];
fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)

scale = prctile(abs(kernels_together(:)),99.97);

subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho.PSV,'relative density PSV','fixed',scale,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.mu.PSV,'relative mu PSV','fixed',scale,stf_PSV);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.lambda.PSV,'relative lambda PSV','fixed',scale,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho.SH,'relative density SH','fixed',scale,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.mu.SH,'relative mu SH','fixed',scale,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative lambda SH','perc',100,stf_PSV);

caxis([-scale scale]);
h = colorbar('horiz');
set(h,'Position',[0.3 0.05 0.4 0.02])