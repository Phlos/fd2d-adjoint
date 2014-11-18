function [diff, fig_diffknl] = diff_kernel_rel(K1, K2)

% plot relative diff of kernels
%
% [diff, fig_diffknl] = diff_kernel_rel(K1, K2)

params = fieldnames(K1);

for i = 1:length(params)
    kaa1 = K1.(params{i}).total;
    kaa2 = K2.(params{i}).total;
%     kaa1(kaa1==0) = eps;
    difftemp = (kaa2 - kaa1);
%     kaa1(kaa1==0) = eps;
    diff.(params{i}).total = difftemp ./ kaa1;
    
%     diffie = diff.(params{i}).total;
%     minnoninf = min(isfinite(diffie));
%     diffie(diffie==-Inf) = minnoninf;
    
    
end

% diff.rho.total = (K2.rho.total - K1.rho.total) ./ K1.rho.total;
% diff.mu.total = (K2.mu.total - K1.mu.total) ./ K1.mu.total;
% diff.lambda.total = (K2.lambda.total - K1.lambda.total) ./ K1.lambda.total;



fig_diffknl = plot_kernels(diff);

end