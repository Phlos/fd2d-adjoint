

function [Params] = update_model(varargin)

% function that calculates the updated model from the (relative!!) kernel,
% the step length and the previous model, or from original input params.
%
% case I:   [Params] = update_model();
%           the func will determine the original parameters from the file
%           input_parameters
% case II:  [Params] = update_model(K_rel, step, Params_previous);
%           the function will calculate a new, smoothed model based on step
%           length, (relative!) kernel and previous model
%
% INPUT (in case II):
% - K_rel:  relative kernel (struct) -- needs K_rel.{parameter}.total at
%           the very least
% - step:   step length
% - Params_previous: the previous model parameter values. Struct with
%   Params_previous.rho, Params_previous.mu and Parama_previous.lambda.
%
% OUTPUT:
% - the new, smoothed updated model. Struct with Params.mu, Params.rho,
%   Params.lambda
%

    Params = checkargs(varargin(:));

end

function Params = checkargs(args)

narg = length(args);


% if no input, get input from file
if narg == 0;
    
    % initialise original parameters
    input_parameters;
    [mu,rho,lambda]=define_material_parameters(nx,nz,model_type);
    Params.rho = rho;
    Params.mu = mu;
    Params.lambda = lambda;
    
    
% otherwise, update model using step length
elseif narg == 3;
    
%     args = varargin(:);
    
    % read input arguments
    K_rel = args{1};
    step = args{2};
    Params_previous = args{3};
 
    % calculate model update
    Params_sharp.rho = Params_previous.rho .* (1 - step * K_rel.rho.total);
    Params_sharp.mu = Params_previous.mu .* (1 - step * K_rel.mu.total);
    Params_sharp.lambda = Params_previous.lambda .* (1 - step * K_rel.lambda.total);
    
    % smooth the updated model!
    Params = smooth_model(Params_sharp, 11, 5);
    
% otherwise error
else
    
    disp(['narg = ',num2str(narg)]);
    error('the number of input variables is not 0 or 3');
    
end

end