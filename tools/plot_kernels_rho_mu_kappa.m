% plot kernels for SH and P-SV in rho, mu, lambda


kernels = figure;
% P-SV
subplot(3,2,1)
plot_kernel(X,Z,K.rho1.PSV,'rho - P-SV',99.95,stf_PSV);
subplot(3,2,3)
plot_kernel(X,Z,K.mu1.PSV,'mu - P-SV',99.95,stf_PSV);
subplot(3,2,5)
plot_kernel(X,Z,K.kappa1.PSV,'kappa - P-SV',99.95,stf_PSV);
% SH
subplot(3,2,2)
plot_kernel(X,Z,K.rho1.SH,'rho - SH',99.97,stf_PSV);
subplot(3,2,4)
plot_kernel(X,Z,K.mu1.SH,'mu - SH',99.98,stf_PSV);
subplot(3,2,6)
plot_kernel(X,Z,zeros(size(X')),'kappa - SH',100,stf_PSV);
