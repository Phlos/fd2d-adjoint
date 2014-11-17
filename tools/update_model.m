

function [Model] = update_model(varargin)

% function that calculates the updated model from the (relative!!) kernel,
% the step length and the previous model, or from original input params.
% Update is carried out in rho-mu-lambda parametrisation!
%
% SYNTAX:
% case I:   [Model] = update_model();
%           the func will determine the original parameters from the file
%           input_parameters
% case II:  [Model] = update_model(K_rel, step, Model_previous);
%           the function will calculate a new, smoothed model based on step
%           length, (relative!) kernel and previous model
%
% INPUT (in case II):
% - K_rel:  relative kernel (struct) -- needs K_rel.{parameter}.total at
%           the very least
% - step:   step length
% - Model_previous: the previous model parameter values. Struct with
%   Model_previous.rho, Model_previous.mu and Model_previous.lambda.
%
% OUTPUT:
% - the new, smoothed updated model. Struct with Model.mu, Model.rho,
%   Model.lambda
%

path(path,'../input')
path(path,'../code')
path(path,'../code/propagation')

    Model = checkargs(varargin(:));

end

function Model = checkargs(args)

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

narg = length(args);


% if no input, get input from file
if narg == 0;
    
    % initialise original parameters
    input_parameters;
    [Model.mu,Model.rho,Model.lambda]=define_material_parameters(nx,nz,model_type);
%     Model.mu = mu;
%     Model.lambda = lambda;
%     Model.rho = rho;
    
% if input is a number and a model number at that:    
elseif narg == 1;
    modelnr = args{1};
    
    [Model.mu,Model.rho,Model.lambda] = ...
                              define_material_parameters(nx,nz,modelnr);
%     Model.rho = rho;
%     Model.mu = mu;
%     Model.lambda = lambda;
    
% otherwise, update model using step length
elseif narg == 3;
    
%     args = varargin(:);
    
    % read input arguments
    K_rel = args{1};
    step = args{2};
    Model_previous = args{3};
    
%     which K_rel
%     K_rel
    whos K_rel.rho.total
    
%     K_sm = smooth_kernels(K_rel, smoothnp, smoothgwid);
    K_sm.rho =      filter_kernels(K_rel.rho.total,smoothgwid);
    K_sm.mu =       filter_kernels(K_rel.mu.total,smoothgwid);
    K_sm.lambda =   filter_kernels(K_rel.lambda.total,smoothgwid);
 
    % calculate model update
    Model.mu =     Model_previous.mu .* (1 - step * K_sm.mu);
    Model.lambda = Model_previous.lambda .* (1 - step * K_sm.lambda);
    Model.rho =    Model_previous.rho .* (1 - step * K_sm.rho);
    
    % smooth the updated model? --> can't do that because then the model

   
% otherwise error
else
    
    disp(['narg = ',num2str(narg)]);
    error('the number of input variables is not 0 or 3');
    
end

end