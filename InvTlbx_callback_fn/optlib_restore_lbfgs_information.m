function [previous_models, previous_gradients] = optlib_restore_lbfgs_information(usr_par)
% OPTLIB_RESTORE_LBFGS_INFORMATION This function allows to warm-start LBFGS
% by reusing the information from previous iterations. If the function
% returns empty structs for previous_models and previous_gradients,
% LBFGS will start with a steepest decent step.
%
% INPUT:
% usr_par : auxiliary user defined parameters (optional)
%
% OUTPUT:
% previous_models : matrix with the models from previous iterations stored
% columnwise with the last column containing the most recent model
%
% previous_gradients : matrix with the gradients from previous iterations stored
% columnwise with the last column containing the most recent gradient. Must
% have the same size as previous_models.
%

% number of models and gradients stored
num_mod_grad_stored = 5;

% use all but the last model and gradient in usr_par:
% (because usr_par.Model(end) = current model
%      and usr_par.K_rel(end) = current gradient )

% determine length of usr_par.Model, usr_par.Kernel
nprev = numel(usr_par.Model) - 1;

% loop over the last 5 models to fill in previous models
if nprev > 0
    maxprev = min(nprev, num_mod_grad_stored);
    for iprev = 1: maxprev
        
        % previous model and gradient under consideration
        Mod_iprev = usr_par.Model(end - iprev);
        Knl_iprev = usr_par.K_abs(end - iprev);
        
        % load model nr countback into m
        m = map_parameters_to_m(Mod_iprev, usr_par);
        
        usr_par.Mod_current = usr_par.Model(end);
        gm = map_gradparameters_to_gradm(Knl_iprev, usr_par);
        
        % add m to previous_models
        previous_models(:,maxprev - iprev + 1) = m;
        previous_gradients(:,maxprev - iprev + 1) = gm;
    end
    
else
    previous_models=[];
    previous_gradients=[];
end
    
end
