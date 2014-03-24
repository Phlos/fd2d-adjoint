% plot relative kernels rho mu kappa

% figure position and size
kernels_together = [K_rel.rho2.PSV; K_rel.rho2.SH; ...
                    K_rel.vs2.PSV; K_rel.vs2.SH; ...
                    K_rel.vp2.PSV];
fig_knl = figure;
set(fig_knl,'OuterPosition',pos_knl)
                
scale = prctile(kernels_together(:),99.97);

% actual plotting
subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho2.PSV,'relative density PSV','fixed',scale,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.vs2.PSV,'relative v_s PSV','fixed',scale,stf_PSV);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.vp2.PSV,'relative v_p PSV','fixed',scale,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho2.SH,'relative density SH','fixed',scale,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.vs2.SH,'relative v_s SH','fixed',scale,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative v_p SH','perc',100,stf_PSV);