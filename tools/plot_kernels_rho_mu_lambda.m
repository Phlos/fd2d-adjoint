% plot kernels for SH and P-SV in rho, mu, lambda

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
