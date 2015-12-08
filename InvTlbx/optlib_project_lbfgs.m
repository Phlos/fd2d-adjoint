
function [mp] = optlib_project_lbfgs(m, A, b, LBFGS_data)

 HAt = optlib_apply_lbfgs_inverse_hessian(A', LBFGS_data);
 AAt_inv = inv(A * HAt);
 mp = m - (HAt * AAt_inv * A) * m + HAt * AAt_inv * b;

end


% 
% function [rho_p] = project(rho, A, b )
%     AAt_inv = inv(A * A');
%     rho_p = rho(:) - (A' * AAt_inv * A) * rho(:) + A * AAt_inv * b;
% end
% 
% 
% 
% function [rho_p] = project(rho, A, b, func )
%     
%     HAt = func(At);
%     AAt_inv = inv(A * HAt);
%     rho_p = rho(:) - (HAt * AAt_inv * A) * rho(:) + HAt * AAt_inv * b;
% end