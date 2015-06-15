 
function [jm] = eval_objective(m, usr_par)
% EVAL_OBJECTIVE function to evaluate the objective function at a given
% model m.
%
% Input:
% m : model
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% jm : objective value (double)
%
% See also EVAL_GRAD_OBJECTIVE and EVAL_OBJECTIVE_AND_GRADIENT.

%% initialise stuff
misfit_init = usr_par.misfit_init;
whichFrq    = usr_par.whichFrq;
g_obs       = usr_par.g_obs;
sEventInfo  = usr_par.sEventInfo;
sEventObs   = usr_par.sEventObs;
InvProps    = usr_par.InvProps;


%% convert variable structures InvTbx -> my stuff
[Model] = map_m_to_parameters(m, usr_par);

%% calculate misfits
disp(['calculating current misfit']);

[misfit_total, misfit_seis, misfit_grav, ...
        giter, g_src, sEventRecIter, sEventAdstfIter] = calc_misfits(Model, ...
                  g_obs, misfit_init(whichFrq).grav , ...
                  sEventInfo, sEventObs, misfit_init(whichFrq).seis, ...
                  'nosavefields','yessaveplots');

% InvProps.misfit(iter) = misfit_total;
% InvProps.misfitseis(iter) = misfit_seis;
% InvProps.misfitgrav(iter) = misfit_grav;

%% OUTPUT to Inversion Toolbox structure
jm = misfit_total;

usr_par.g_src           = g_src;
usr_par.sEventRecIter   = sEventRecIter;
usr_par.sEventAdstfIter = sEventAdstfIter;

% usr_par.InvProps = InvProps;
% usr_par.misfit_total  = misfit_total;
% usr_par.misfit_seis   = misfit_seis;
% usr_par.misfit_grav   = misfit_grav;

end