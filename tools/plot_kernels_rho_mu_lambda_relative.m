% plot relative kernels rho mu lambda

figure;
subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho.PSV,'relative density PSV',99.97,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.mu.PSV,'relative mu PSV',99.97,[1 0]);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.lambda.PSV,'relative lambda PSV',99.97,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho.SH,'relative density SH',99.97,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.mu.SH,'relative mu SH',99.98,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative lambda SH',99.97,stf_PSV);