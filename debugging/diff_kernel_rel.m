function [diff, fig_diffknl] = diff_kernel_rel(K1, K2)

% plot relative diff of kernels
%
% [diff, fig_diffknl] = diff_kernel_rel(K1, K2)

params = fieldnames(K1);

for i = 1:length(params)
    diff.(params{i}).total = (K2.(params{i}).total - K1.(params{i}).total) ./ K1.(params{i}).total;
end

% diff.rho.total = (K2.rho.total - K1.rho.total) ./ K1.rho.total;
% diff.mu.total = (K2.mu.total - K1.mu.total) ./ K1.mu.total;
% diff.lambda.total = (K2.lambda.total - K1.lambda.total) ./ K1.lambda.total;
fig_diffknl = plot_kernels(diff);

end