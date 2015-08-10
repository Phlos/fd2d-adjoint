function K_rel = calculate_relative_kernels(K, varargin)

% calculates kernels relative to a model (K_rel = K * modelparameter)
%
% SYNTAX:
% - K_rel = calculate_relative_kernels(K)
% - K_rel = calculate_relative_kernels(K, modelnr)
% - K_rel = calculate_relative_kernels(K, model)
%
% INPUT:
% - K       the kernels, struct of shape K.{parameter}.{component}
% - modelnr the model number of one of the standard models -- this model
%           will then be calculated with the help of input_parameters
% - Model   the model, struct of shape Model.rho, .mu, .lambda
% --> if no modelnr of Model is supplied, the modelnr will be determined
%     from the starting model defined in input_parameters.
%
% OUTPUT:
% - K_rel   the kernels of the same shape as K, but now relative to the
%           model. [K_rel = K * modelparameter]

% get the model in all parametrisations
Model = checkargs(varargin(:));
if isfield(Model, 'rho') && isfield(Model, 'vs') && isfield(Model, 'vp')
    Model.rho2 = Model.rho;
    Model.vs2 = Model.vs;
    Model.vp2 = Model.vp;
elseif isfield(Model, 'rho') && isfield(Model, 'mu') && isfield(Model, 'lambda')
    Model_rvv = change_parametrisation('rhomulambda','rhovsvp',Model);
    Model.rho2 = Model_rvv.rho;
    Model.vs2 = Model_rvv.vs;
    Model.vp2 = Model_rvv.vp;
end
        

% loop over model parameters
fn_params = fieldnames(K);
for i = 1:length(fn_params)

    % loop over parameter components
    fn_comps = fieldnames(K.(fn_params{i}));
    for j = 1:length(fn_comps)
        
        % calculate the relative kernel
        K_rel.(fn_params{i}).(fn_comps{j}) = K.(fn_params{i}).(fn_comps{j}) .* Model.(fn_params{i});

    end
end


end

function Model = checkargs(args)

% either the model from input_parameters or a model number or a supplied
% model
%
% INPUT
% - option 1: nothing --> get from input_parameters
% - option 2: model number --> calculate model based on model number
% - option 3: input model --> pass on

% size(args)
narg = length(args);


% initialise domain
input_parameters;

% if no input, get from input_parameters
if narg == 0;
%     disp 'we detected no input... getting model from input_parameters'
    [Model.mu,Model.rho, Model.lambda]=define_material_parameters(nx,nz,model_type);
    
    
 
elseif narg == 1;
%     disp 'one argument, it seems'

    % if input is a number and a model number at that:   
    if isfloat(args{1})
%         disp 'we detected a model number'
        modelnr = args{1};
        [Model.mu,Model.rho,Model.lambda] = ...
                              define_material_parameters(nx,nz,modelnr);
%         Model
                          
    % if input is a struct it is the struct with a model in it
    elseif isstruct(args{1})
%         disp 'we detected a structure array'
        Model = args{1};

    else
        banaan = args{1};
        which('banaan')
        whos('banaan')
        error('the input variable type was not recognised')
    end
else
    error('the input to calculate_relative_kernels seems flawed');
end
end