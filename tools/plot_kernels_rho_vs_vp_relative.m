% plot relative kernels rho mu kappa

figure;
subplot (3,2,1);
plot_kernel(X,Z,K_rel.rho2.PSV,'relative density PSV',99.97,stf_PSV);
subplot (3,2,3);
plot_kernel(X,Z,K_rel.vs2.PSV,'relative v_s PSV',99.97,[1 0]);
subplot (3,2,5);
plot_kernel(X,Z,K_rel.vp2.PSV,'relative v_p PSV',99.97,stf_PSV);

subplot (3,2,2);
plot_kernel(X,Z,K_rel.rho2.SH,'relative density SH',99.97,stf_PSV);
subplot (3,2,4);
plot_kernel(X,Z,K_rel.vs2.SH,'relative v_s SH',99.98,stf_PSV);
subplot (3,2,6);
plot_kernel(X,Z,zeros(size(X')),'relative v_p SH',99.97,stf_PSV);