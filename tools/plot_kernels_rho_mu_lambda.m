% plot kernels for SH and P-SV in rho, mu, lambda


kernels = figure;
% P-SV
subplot(3,2,1)
plot_kernel(X,Z,K.rho.PSV,'rho - P-SV',99.95,stf_PSV);
subplot(3,2,3)
plot_kernel(X,Z,K.mu.PSV,'mu - P-SV',99.95,stf_PSV);
subplot(3,2,5)
plot_kernel(X,Z,K.lambda.PSV,'lambda - P-SV',99.95,stf_PSV);
% SH
subplot(3,2,2)
plot_kernel(X,Z,K.rho.SH,'rho - SH',99.97,stf_PSV);
subplot(3,2,4)
plot_kernel(X,Z,K.mu.SH,'mu - SH',99.995,stf_PSV);
subplot(3,2,6)
plot_kernel(X,Z,zeros(size(X')),'lambda - SH',100,stf_PSV);
