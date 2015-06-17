function [Kernel] = map_gradm_to_gradparameters(gm, usr_par)
% MAP_GRADM_TO_GRADPARAMETERS function to map the model gradient gm to 
% physical parameter gradients rho, lambda and mu. (RELATIVE GRADIENT)
%
% Input: 
% gm
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% Kernel
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.

%% input 
% input
input_parameters;

parametrisation = usr_par.parametrisation;

%% cut gradient 

% determine which part of the m vector is which parameters
switch parametrisation
    case 'rhomulambda'
        
        % cut gradient into 3 params
        gm1 = gm(          1 :   nx*nz);
        gm2 = gm(  nx*nz + 1 : 2*nx*nz);
        gm3 = gm(2*nx*nz + 1 : 3*nx*nz);
        
        % reshape
        Kernel.rho.total     = reshape(gm1, nx, nz);
        Kernel.mu.total      = reshape(gm2, nx, nz);
        Kernel.lambda.total  = reshape(gm3, nx, nz);
        
    case 'rhovsvp'
        
        gm4 = m(          1 :   nx*nz);
        gm5 = m(  nx*nz + 1 : 2*nx*nz);
        gm6 = m(2*nx*nz + 1 : 3*nx*nz);
        
        % reshape
        Kernel.rho2.total = reshape(gm4, nx, nz);
        Kernel.vs2.total  = reshape(gm5, nx, nz);
        Kernel.vp2.total  = reshape(gm6, nx, nz);
        
end





end