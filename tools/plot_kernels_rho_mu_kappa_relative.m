% plot relative kernels rho mu kappa

figure;
subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho1.PSV,'relative density PSV',99.97,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.mu1.PSV,'relative mu PSV',99.97,[1 0]);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.kappa1.PSV,'relative kappa PSV',99.97,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho1.SH,'relative density SH',99.97,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.mu1.SH,'relative mu SH',99.98,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative kappa SH',99.97,stf_PSV);