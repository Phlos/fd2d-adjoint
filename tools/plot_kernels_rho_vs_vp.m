% plot kernels for SH and P-SV in rho, mu, lambda


kernels = figure;
% P-SV
subplot(3,2,1)
plot_kernel(X,Z,K.rho2.PSV,'rho - P-SV',99.95,stf_PSV);
subplot(3,2,3)
plot_kernel(X,Z,K.vs2.PSV,'vs - P-SV',99.95,stf_PSV);
subplot(3,2,5)
plot_kernel(X,Z,K.vp2.PSV,'vp - P-SV',99.95,stf_PSV);
% SH
subplot(3,2,2)
plot_kernel(X,Z,K.rho2.SH,'rho - SH',99.97,stf_PSV);
subplot(3,2,4)
plot_kernel(X,Z,K.vs2.SH,'v_s - SH',99.98,stf_PSV);
subplot(3,2,6)
plot_kernel(X,Z,zeros(size(X')),'v_p - SH',100,stf_PSV);
