function [K_abs] = map_gradm_to_gradparameters(gm, usr_par)
% MAP_GRADM_TO_GRADPARAMETERS function to map the model gradient gm to 
% physical parameter gradients rho, lambda and mu. (RELATIVE GRADIENT)
%
% Input: 
% gm
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% Kernel : the absolute kernels
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.

%% input 
% input
input_parameters;

parametrisation = usr_par.parametrisation;

%% reparametrise gm to relative kernels in physical parameters

% fixing vs and vp: there are only 1*nx*nz free parameters: rho
if strcmp(fix_velocities,'yes')
    
    if strcmp(parametrisation, 'rhomulambda')
       
        K_rel.rho.total     = reshape(gm, nx, nz);
        K_rel.mu.total      = zeros(size(K_rel.rho.total));
        K_rel.lambda.total  = zeros(size(K_rel.rho.total));
    else
	    K_rel.rho2.total = reshape(gm, nx, nz);
	    K_rel.vs2.total = zeros(size(K_rel.rho2.total));
	    K_rel.vp2.total = zeros(size(K_rel.rho2.total));

%    	    error('parametrisation must be rhomulambda if fixing velocities');
    end

% fixing density: there are 2*nx*nz free parameters: mu/lambda or vs/vp
elseif strcmp(fix_density, 'yes')
	
	switch parametrisation

	case 'rhomulambda'

	    gm2 = gm(          1 :   nx*nz);
            gm3 = gm(  nx*nz + 1 : 2*nx*nz);

	    % reshape
		K_rel.mu.total     = reshape(gm2, nx, nz);
		K_rel.lambda.total = reshape(gm3, nx, nz);
		K_rel.rho.total    = zeros(size(K_rel.mu.total));

	case 'rhovsvp'

	    gm5 = gm(          1 :   nx*nz);
            gm6 = gm(  nx*nz + 1 : 2*nx*nz);
            
            % reshape
            K_rel.vs2.total  = reshape(gm5, nx, nz);
            K_rel.vp2.total  = reshape(gm6, nx, nz);
            K_rel.rho2.total = zeros(size(K_rel.vs2.total));

    	end

% no fixing of parameters: there are 3*nx*nz free parameters
else

    % determine which part of the m vector is which parameters
    switch parametrisation
        case 'rhomulambda'
            
            % cut gradient into 3 params
            gm1 = gm(          1 :   nx*nz);
            gm2 = gm(  nx*nz + 1 : 2*nx*nz);
            gm3 = gm(2*nx*nz + 1 : 3*nx*nz);
            
            % reshape
            K_rel.rho.total     = reshape(gm1, nx, nz);
            K_rel.mu.total      = reshape(gm2, nx, nz);
            K_rel.lambda.total  = reshape(gm3, nx, nz);
            
        case 'rhovsvp'
            
            gm4 = gm(          1 :   nx*nz);
            gm5 = gm(  nx*nz + 1 : 2*nx*nz);
            gm6 = gm(2*nx*nz + 1 : 3*nx*nz);
            
            % reshape
            K_rel.rho2.total = reshape(gm4, nx, nz);
            K_rel.vs2.total  = reshape(gm5, nx, nz);
            K_rel.vp2.total  = reshape(gm6, nx, nz);
            
    end

end

K_abs = calculate_unrelative_kernels(K_rel, usr_par.Model_bg);
K_abs = change_parametrisation_kernels(parametrisation, 'rhomulambda', K_abs, usr_par.Model_bg);


end
