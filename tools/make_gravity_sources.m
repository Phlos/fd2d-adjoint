function [g_src, misfit_g] = make_gravity_sources(g, g_obs, which_grav, varargin)

% calculates the gravity sources & gravity misfit based on current &
% observed gravity field. 
% - two-norm of (g_rec - g_obs) (full gravity vector)
% - two-norm of (W_rec - W_obs) (gravity potential)
%
% [g_src, misfit_g] = make_gravity_sources(g, g_obs, which_grav, varargin)
%
% INPUT:
% - g: struct with .x, .z, .mag, .pot
%      where x and z form the gravity vector, mag their magnitude and .pot
%      is the gravity potential at the receiver location
% - g_obs: same, but the observed values
% - which_grav: the switch which determines which datum is used as the
%   actual gravity misfit source
%   * two-norm of full gravity vector, or
%   * 'two-norm' of scalar gravity potential
%
% OUTPUT:
% - g_src:  the 'gravity source' that is used in calculating the gravity
%           kernels later (in compute_kernels_gravity)
%           .x, .z and .pot are calculated
% - misfit_g: the gravity misfit, struct with .x .z .pot
%           used to calculate the full misfit (duh.)
%
% NOTE:
% I could also construct a gravity misfit based on the _magnitude_ of the
% gravity vector g.mag. just add two lines of code for calculating
% g_src.mag and misfit_g.mag...
%
% -- Nienke Blom, 24-9-2015

% obtain scalingfactor from varargin
% scalingfactor = checkargs(varargin(:));

%% misfit  |g_rec - g_obs|^2

% gravity source (g vector)
g_src.x = g.x - g_obs.x;
g_src.z = g.z - g_obs.z;

% quadratic misfit (g vector)
misfit_g.x = sum(g_src.x .* g_src.x);
misfit_g.z = sum(g_src.z .* g_src.z);


%% misfit (gpotential_rec - gpotential_obs)^2
% (W_rec - W_obs)^2

% gravity source (g potential)
g_src.pot = g.pot - g_obs.pot;
% quadratic misfit (g potential)
misfit_g.pot = sum( g_src.pot .* g_src.pot );


%% determine actual, total misfit
% which_grav is used here to determine misfit.total

% if(exist('which_grav', 'var')
if strcmp(which_grav, 'g_potential')
    misfit_g.total = misfit_g.pot;
elseif strcmp(which_grav, 'g_vector')
    misfit_g.total = misfit_g.x + misfit_g.z;
else
    error('which_grav misfit type not recognised');
end
% else
%     warning('ATTENTION!! gravity misfit being used is g_vector');
%     misfit_g.total = misfit_g.x + misfit_g.z;
% end

misfit_g.normd = NaN;
% determine scalingfactor if the misfit has to be scaled by itself
% (iter 1)
% if isnan(scalingfactor)
%     scalingfactor = misfit_g.total;
% end

% % divide-by-zero catch
% if scalingfactor == 0
%     misfit_g.normd = misfit_g.total;
% else
%     misfit_g.normd = misfit_g.total / scalingfactor;
% end


end

% function sf = checkargs(args)
% 
% % checking input to make_gravity_sources
% % no arguments:     sf = 1 (no scaling)
% % args{1} = NaN:    sf = NaN, so that the misfit will be scaled by itself
% % args{1} = number: sf = args{1}, misfit will be scaled by that number.
% 
% nargs = length(args);
% 
% if (nargs == 0)
%     disp 'make_gravity_sources: no scaling factor supplied'
%     sf = 1;
% elseif (nargs == 1 && isnumeric(args{1}))
%     sf = args{1};
% %     disp(['make_gravity_sources: scaling factor = ', num2str(sf)]);
%     if isnan(sf)
% %         disp ' --> misfit will be scaled by itself'
%     end
% else
%     error('too many or wrong input arguments to make_gravity_sources')
% end
% 
% end