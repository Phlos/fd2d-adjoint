function [g_src, misfit_g] = make_gravity_sources(g, g_obs)

%% misfit  |g_rec - g_obs|^2

% gravity source --> this could be extended if we want a different misfit functional
g_src.x = g.x - g_obs.x;
g_src.z = g.z - g_obs.z;

% quadratic misfit.
misfit_g.x = sum(g_src.x .* g_src.x);
misfit_g.z = sum(g_src.z .* g_src.z);
misfit_g.total = misfit_g.x + misfit_g.z;



end