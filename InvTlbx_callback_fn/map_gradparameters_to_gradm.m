
function [gradm] = map_gradparameters_to_gradm(gradrho, gradlambda, gradmu)
% MAP_GRADPARAMETERS_TO_GRADM function to transform the partial derivatives
% w.r.t. the physical parameters rho, lambda to the partial derivatives
% w.r.t. to the model variable m. This requires to apply the chain rule
% which is implicitly given by the function MAP_M_TO_PARAMETERS.
%
% Input: 
% gradrho 
% gradlambda
% gradmu
%
% Output:
% gradm
%
% See also MAP_M_TO_PARAMETERS.

end