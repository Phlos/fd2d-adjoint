function [rho_p, fig_rhop] = project_constraints_sd(X, Z, rho, rec_g, g_obs, props_obs, which_constraints)

% function project_constraints_sd
%
% input:
% - rho: density field
% - X:   X location field
% - Z:   Z location field
% - g_obs: observed gravity values
% - props_obs: observed mass/moi properties of model
% - rec_g: locations of gravity receivers
% - which_constraints, a cell array which may include:   
%                                   'g_vec', 'g_pot'
%
% output:
% - rho_p:  projected density field
% - grav_proj: gravity field of projected density field
% - grav_diff: normalised difference between input & output gravity fields
% - fig_rhop:  figure with the original, projected and diff rho fields
%

%% preparation

%  0a) test if all fields are in the right shape
if ~( isequal(size(rho), size(X)) && isequal(size(rho), size(Z)) )
    Xn = X'; Zn = Z';
    clearvars X Z;
    if ~( isequal(size(rho), size(Xn)) && isequal(size(rho), size(Zn)) )
        error('X, Z, rho fields are not of the same size');
    end
else
    Xn = X; Zn = Z;
end

% make column vector out of rho
rho_vec = rho(:);

%% initialise constraints

 [c_A, c_b] = init_constraints(Xn, Zn, rec_g, g_obs, props_obs, which_constraints);

%% projection
% rho_p = rho - A'* (A A')^-1 * (A * rho - b)
AAA = c_A' * (c_A*c_A')^-1;
Ar_u_minb = c_A * rho_vec - c_b;
update = AAA * Ar_u_minb;
rho_p_vec = rho_vec - update;


%% output
% props_proj = NaN;
% props_diff = NaN;
%- rho in reshaped form
nx = size(rho, 1); nz = size(rho,2);
rho_p = reshape(rho_p_vec, nx, nz);

% gravity field of output rho_p (to check)
g_proj = calculate_gravity_field(rho_p, rec_g, 'noplot');
gravs = fieldnames(g_obs);
% for ii = 1:numel(gravs)
%     g_diff.(gravs{ii}) = (g_proj.(gravs{ii}) - g_obs.(gravs{ii})) ./ g_obs.(gravs{ii});
% end

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
title('original density field');
subplot(3,1,2); colorbar;
title('new density field');
subplot(3,1,3); colorbar;
title('difference density field');

end