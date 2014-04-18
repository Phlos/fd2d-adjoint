

function [Params] = update_model(K,step)

% initialise original parameters
input_parameters;
[mu,rho,lambda]=define_material_parameters(nx,nz,model_type);


% update model using step length

Params.rho = rho + step * K.rho.total;
Params.mu = mu + step * K.mu.total;
Params.lambda = lambda + step * K.lambda.total;

end