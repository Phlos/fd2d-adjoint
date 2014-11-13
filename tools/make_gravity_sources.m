function [g_src, misfit_g] = make_gravity_sources(g, g_obs)

% gravity source --> this could be extended if we want a different misfit functional
g_src.x = g.x - g_obs.x;
g_src.z = g.z - g_obs.z;

% quadratic misfit.
misfit_g = sum(g_src.x .* g_src.x  +  g_src.z .* g_src.z);

end