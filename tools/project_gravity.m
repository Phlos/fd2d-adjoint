function [rho_p, g_proj, g_diff, fig_rhop] = project_gravity(rho, X, Z, g_obs, rec_g, options)

% function project_gravity
%
% input:
% - rho: density field
% - X:   X location field
% - Z:   Z location field
% - g_obs: observed gravity values
% - options:    .which_constraints, a cell array which may include:   
%                                   'g_vec', 'g_pot'
%
% output:
% - rho_p:  projected density field
% - grav_proj: gravity field of projected density field
% - grav_diff: normalised difference between input & output gravity fields
% - fig_rhop:  figure with the original, projected and diff rho fields
%

%- 0) preparation:
nrec = numel(rec_g.x);

%  0a) test if all fields are in the right shape
Xn = X'; Zn = Z';
clearvars X Z;
if ~( isequal(size(rho), size(Xn)) && isequal(size(rho), size(Zn)) )
    error('X, Z, rho fields are not of the same size');
end
if ~(isequal(nrec, numel(g_obs.mag)) )
    error('g_obs and g_rec have different/wrong sizes')
end

%  0b) determine constants G, dx, dz
G = 6.67384e-11;        % universal gravity const. 6.67384 * 10-11 m3 kg-1 s-2
dx = Xn(2,1) - Xn(1,1);
dz = Zn(1,2) - Zn(1,1);
% NOTE: dx and dz are strictly not necessary in these calculations, so I
%       might as well leave them out. It just feels weird not to use them.


%  0c) determine which constraints
%      'g_vec' or 'g_pot'


%  0d) determine how many constraints
% if strcmp(constr_type, 'g_vec')
%     nconstr = 2*nrec
% elseif strcmp(constr_type, 'g_pot')
    % gravity potential
    nconstr = nrec;
% end

%- 1) transform rho to 1D vector
rho_vec = rho(:);

for ii = 1:nrec
    %- 2a) calculate the gravity influence matrices per receiver
    rX = rec_g.x(ii) - Xn;
    rZ = rec_g.z(ii) - Zn;
    rLen = sqrt(rX .^2 + rZ .^2);
    
    % gravity vector
    ArX = -G ./ rLen .^3  .* rX;
    ArZ = -G ./ rLen .^3  .* rZ;
    % gravity potential
    ArP  = -G ./ rLen;
    
    %- 2b) transform the gravity influence matrix components to 1D vectors
    ArX_vec = ArX(:);
    ArZ_vec = ArZ(:);
    ArP_vec = ArP(:);

    %- 3a) build up A matrix
    % gravity potential
    A(ii,:) = ArP_vec * dx * dz;
    
    %- 3b) build up b constraints vector;
    % gravity potential
    b(ii) = g_obs.pot(ii);
end



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

% gravity field of output rho_p (to check)
g_proj = calculate_gravity_field(rho_p, rec_g, 'noplot');
gravs = fieldnames(g_obs);
for ii = 1:numel(gravs)
    g_diff.(gravs{ii}) = (g_proj.(gravs{ii}) - g_obs.(gravs{ii})) ./ g_obs.(gravs{ii});
end

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