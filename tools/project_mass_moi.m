function [rho_p, props_proj, props_diff, fig_rhop] = project_mass_moi(rho, X, Z, constraints, options);

% function project_mass_moi
%
% input:
% - rho: density field
% - X:   X location field
% - Z:   Z location field
% - constraints:    .mass: the total mass of the model
%                   .I_11: the (1,1) component of the inertia tensor
%                   .I_12: the (1,2) component of the inertia tensor
%                   .I_22: the (2,2) component of the inertia tensor
% - options:    .which_constraints, a cell array which may include:   
%                                   'mass', 'I_11', 'I_12', 'I_22'
%
% output:
% - rho_p:  projected density field
%

%- 0) preparation:

%  0a) test if all fields are in the right shape
Xn = X'; Zn = Z';
clearvars X Z;
% if (size(rho) ~= size(X)) || (size(rho) ~= size(Z))
if ~( isequal(size(rho), size(Xn)) && isequal(size(rho), size(Zn)) )
    error('fields are not of the same size');
end

%  0b) determine dx, dz
dx = Xn(2,1) - Xn(1,1);
dz = Zn(1,2) - Zn(1,1);
% NOTE: dx and dz are strictly not necessary in these calculations, so I
%       might as well leave them out. It just feels weird not to use them.

%  0c) determine which constraints


%- 1) transform X, Z, rho to 1D vectors
rho_vec = rho(:);

%- 2a) calculate the inertia matrix R_ij
R_11 = Zn .^2;
R_12 = -Xn .* Zn;
R_22 = Xn .^2;

%- 2b) transform the inertia matrix components to 1D vectors
R_11_vec = R_11(:); % transpose is to make it a row vector
R_12_vec = R_12(:);
R_22_vec = R_22(:);

%- 3a) build up A matrix
A(1,:) = ones(size(rho_vec')) * dx * dz;
A(2,:) = R_11_vec * dx * dz;
A(3,:) = R_22_vec * dx * dz;
A(4,:) = R_12_vec * dx * dz;

%- 3b) build up b constraints vector;
b(1) = constraints.mass;
b(2) = constraints.I_11;
b(3) = constraints.I_22;
b(4) = constraints.I_12;
b = b'; % make column vector

%- 4) perform least squares projection
% rho_p = rho - A'* (A A')^-1 * (A * rho - b)
AAA = A' * (A*A')^-1;
Ar_u_minb = A * rho_vec - b;
update = AAA * Ar_u_minb;
rho_p_vec = rho_vec - update;


%% output
%- rho in reshaped form
nx = size(rho, 1); nz = size(rho,2);
rho_p = reshape(rho_p_vec, nx, nz);

% properties of output rho_p (to check)
props_proj = calc_mass_moi(rho_p, X, Z);
props_diff.mass = (constraints.mass - props_proj.mass) / constraints.mass;
props_diff.I_11 = (constraints.I_11 - props_proj.I_11) / constraints.I_11;
props_diff.I_12 = (constraints.I_12 - props_proj.I_12) / constraints.I_12;
props_diff.I_22 = (constraints.I_22 - props_proj.I_22) / constraints.I_22;

% plot output:
fig_rhop = figure;
subplot(3,1,1); pcolor(Xn,Zn,rho); axis image; shading interp;
subplot(3,1,2); pcolor(Xn,Zn,rho_p); axis image; shading interp;
subplot(3,1,3); pcolor(Xn,Zn,rho_p-rho); axis image; shading interp;
minmaks = [2.6000    5.5665] * 1e3;
% fig_rhop.Children(1).CLim = minmaks;
fig_rhop.Children(2).CLim = minmaks;
fig_rhop.Children(3).CLim = minmaks;
subplot(3,1,1); colorbar;
subplot(3,1,2); colorbar;
subplot(3,1,3); colorbar;

end
