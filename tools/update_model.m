

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
% - K_abs:  relative kernel (struct) -- needs K_abs.{parameter}.total at
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
path(path,'../tools');

    Model = checkargs(varargin(:));
%     Model

end

function Modelout = checkargs(args)

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

narg = length(args);


% if no input, get input from file
if narg == 0;
    
    % initialise original parameters
    input_parameters;
    [Modelo.mu,Modelo.rho,Modelo.lambda] = ...
                              define_material_parameters(nx,nz,model_type);
    
    % this is to make sure that the params are stored in the right order
    Modelout.rho = Modelo.rho;
    Modelout.mu = Modelo.mu;
    Modelout.lambda = Modelo.lambda;
    
% if input is a number and a model number at that:    
elseif narg == 1;
    modelnr = args{1};
    
    [Modelo.mu,Modelo.rho,Modelo.lambda] = ...
                              define_material_parameters(nx,nz,modelnr);
                          
    % this is to make sure that the params are stored in the right order
    Modelout.rho = Modelo.rho;
    Modelout.mu = Modelo.mu;
    Modelout.lambda = Modelo.lambda;
    
% otherwise, update model using step length
elseif narg == 3 || narg == 4
    
%     disp '3 or 4 input to update_model'
        
    % read input arguments
    K_in = args{1};
    step = args{2};
    Model_prev = args{3};
    if narg == 4
        parametrisation = args{4};
    else
        parametrisation = 'rhomulambda';
    end
    
%     disp(['updating model in ', parametrisation])
    
%     disp 'after reading input'
%     parametrisation
%     K_in
%     Model_prev
    
    if strcmp(parametrisation,'rhomulambda')
        K_abs = K_in;
    elseif strcmp(parametrisation,'rhovsvp')
        K_abs = change_parametrisation_kernels('rhomulambda','rhovsvp',K_in, Model_prev);
        Model_prev = change_parametrisation('rhomulambda', 'rhovsvp',Model_prev);
    else
        error('Parametrisation not recognised');
    end
    
%     disp 'after changing the parametrisation of input K and Model';
%     K_abs
%     Model_prev
    
    % calculate relative kernels
    K_rel = calculate_relative_kernels(K_abs, Model_prev);
    
%     disp 'after calculating relative kernels'
%     K_rel
    
    % smooth the kernels
    % K_sm = smooth_kernels(K_rel, smoothnp, smoothgwid);
    if strcmp(parametrisation,'rhomulambda')
        K_sm.rho =      filter_kernels(K_rel.rho.total,smoothgwid);
        K_sm.mu =       filter_kernels(K_rel.mu.total,smoothgwid);
        K_sm.lambda =   filter_kernels(K_rel.lambda.total,smoothgwid);
    elseif strcmp(parametrisation,'rhovsvp')
        K_sm.rho2 =      filter_kernels(K_rel.rho2.total,smoothgwid);
        K_sm.vs2 =       filter_kernels(K_rel.vs2.total,smoothgwid);
        K_sm.vp2 =   filter_kernels(K_rel.vp2.total,smoothgwid);
    else
        error('Parametrisation not recognised (smoothing kernels)');
    end
    
%     disp 'after smoothing the kernels'
%     K_sm
    
    % calculate model update
    if strcmp(parametrisation,'rhomulambda')
        Model.rho =    Model_prev.rho .* (1 - step * K_sm.rho);
        Model.mu =     Model_prev.mu .* (1 - step * K_sm.mu);
        Model.lambda = Model_prev.lambda .* (1 - step * K_sm.lambda);
    elseif strcmp(parametrisation,'rhovsvp')
        Model.rho =    Model_prev.rho .* (1 - step * K_sm.rho2);
        Model.vs =     Model_prev.vs .* (1 - step * K_sm.vs2);
        Model.vp =     Model_prev.vp .* (1 - step * K_sm.vp2);
    else
        error('Parametrisation not recognised (updating model)');
    end
    
%     disp 'after calculating the model update'
%     Model
    
    % if update param != rhomulambda, reparametrise model to rhomulambda
    if strcmp(parametrisation,'rhomulambda')
        Modelout = Model;
    elseif strcmp(parametrisation,'rhovsvp')
        Modelout = change_parametrisation('rhovsvp','rhomulambda',Model);
    else
        error('Parametrisation not recognised (reparametrising to rhomulambda)');
    end
    
%     disp 'after changing it to rho-mu-lambda'
%     Modelout
    

% elseif narg==4;
%     
%     % read input arguments
%     K_abs = args{1};
%     step = args{2};
%     Model_prev = args{3};
%     parametrisation = args{4};
%     
%     % calculate relative kernels
%     K_rel = calculate_relative_kernels(K_abs, Model_prev);
%     
%     % smooth the kernels
%     K_sm.rho =      filter_kernels(K_rel.rho.total,smoothgwid);
%     K_sm.mu =       filter_kernels(K_rel.mu.total,smoothgwid);
%     K_sm.lambda =   filter_kernels(K_rel.lambda.total,smoothgwid);
%     
%     % calculate model update
%     Model.mu =     Model_prev.mu .* (1 - step * K_sm.mu);
%     Model.lambda = Model_prev.lambda .* (1 - step * K_sm.lambda);
%     Model.rho =    Model_prev.rho .* (1 - step * K_sm.rho);
%     
%     
%     error('still add the parametrisation thingie!')
   
% otherwise error
else
    
    disp(['narg = ',num2str(narg)]);
    error('the number of input variables is not 0 or 3');
    
end

end