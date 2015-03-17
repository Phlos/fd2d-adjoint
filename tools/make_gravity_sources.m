function [g_src, misfit_g] = make_gravity_sources(g, g_obs, varargin)

% calculates the gravity sources & gravity misfit based on current &
% observed gravity field. 
% - two-norm of (g_rec - g_obs)  --> currently the only option

% obtain scalingfactor from varargin
scalingfactor = checkargs(varargin(:));

%% misfit  |g_rec - g_obs|^2

% gravity source --> this could be extended if we want a different misfit functional
g_src.x = g.x - g_obs.x;
g_src.z = g.z - g_obs.z;

% quadratic misfit.
misfit_g.x = sum(g_src.x .* g_src.x);
misfit_g.z = sum(g_src.z .* g_src.z);
misfit_g.total = misfit_g.x + misfit_g.z;
% determine scalingfactor if the misfit has to be scaled by itself
% (iter 1)
if isnan(scalingfactor)
    scalingfactor = misfit_g.total;
end

% divide-by-zero catch
if scalingfactor == 0
    misfit_g.normd = misfit_g.total;
else
    misfit_g.normd = misfit_g.total / scalingfactor;
end


end

function sf = checkargs(args)

% checking input to make_gravity_sources
% no arguments:     sf = 1 (no scaling)
% args{1} = NaN:    sf = NaN, so that the misfit will be scaled by itself
% args{1} = number: sf = args{1}, misfit will be scaled by that number.

nargs = length(args);

if (nargs == 0)
    disp 'make_gravity_sources: no scaling factor supplied'
    sf = 1;
elseif (nargs == 1 && isnumeric(args{1}))
    sf = args{1};
%     disp(['make_gravity_sources: scaling factor = ', num2str(sf)]);
    if isnan(sf)
%         disp ' --> misfit will be scaled by itself'
    end
else
    error('too many or wrong input arguments to make_gravity_sources')
end

end