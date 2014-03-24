% plot relative kernels rho mu kappa

% figure position and size
kernels_together = [K_rel.rho1.PSV; K_rel.rho1.SH; ...
                    K_rel.mu.PSV; K_rel.mu1.SH; ...
                    K_rel.kappa1.PSV];
fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)
                
scale = prctile(kernels_together(:),99.97);

subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho1.PSV,'relative density PSV','fixed',scale,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.mu1.PSV,'relative mu PSV','fixed',scale,stf_PSV);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.kappa1.PSV,'relative kappa PSV','fixed',scale,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho1.SH,'relative density SH','fixed',scale,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.mu1.SH,'relative mu SH','fixed',scale,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative kappa SH','perc',100,stf_PSV);