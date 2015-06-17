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

    previous_models=[];
    previous_gradients=[];

end
