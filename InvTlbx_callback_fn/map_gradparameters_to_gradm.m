
function [gradm] = map_gradparameters_to_gradm(K_rel, usr_par)
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

% input
input_parameters;

% fixing vs and vp: there are only 1*nx*nz free parameters: rho
if strcmp(fix_velocities,'yes')
    
    if strcmp(parametrisation, 'rhomulambda')
        gradm = [K_rel.rho.total(:) ...
            + K_rel.mu.total(:) ...
            + K_rel.lambda.total(:) ];
    else
        error('parametrisation must be rhomulambda if fixing velocities');
    end
    
% no fixing of parameters: there are 3*nx*nz free parameters
else

    switch parametrisation
        case 'rhomulambda'
            gradm = [K_rel.rho.total(:); ...
                K_rel.mu.total(:); ...
                K_rel.lambda.total(:)];
        case 'rhovsvp'
            if strcmp(fix_velocities,'yes')
                error('WARNING! no fix vp vs defined here!')
            end
            gradm = [K_rel.rho2.total(:); ...
                K_rel.vs2.total(:); ...
                K_rel.vp2.total(:)];
    end


end

end