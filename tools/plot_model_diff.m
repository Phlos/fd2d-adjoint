function [fig_mod, Model_diff] = plot_model_diff(Model1, Model2, varargin)

% function that plots the difference model between Model 1 and Model 2
% Model_diff = Model1 - Model2

% varargin{:}
outparam = checkargs(varargin);

switch outparam
    case 'rhomulambda'
        Model_diff.mu = Model1.mu - Model2.mu;
        Model_diff.lambda = Model1.lambda - Model2.lambda;
        Model_diff.rho = Model1.rho - Model2.rho;
    case 'rhovsvp'
        Model1 = change_parametrisation('rhomulambda','rhovsvp',Model1);
        Model2 = change_parametrisation('rhomulambda','rhovsvp',Model2);
        Model_diff.vs = Model1.vs - Model2.vs;
        Model_diff.vp = Model1.vp - Model2.vp;
        Model_diff.rho = Model1.rho - Model2.rho;
    otherwise
        error('parametrisation not recognised');
end


% plotting the model with the '0' meaning that the white of the plot is
% the value zero.
fig_mod = plot_model(Model_diff,0,varargin{:});


end




function [outparam] = checkargs(arg)

% defines the output parametrisation based on whether any is given (default
% is 'rhomulambda')

narg = length(arg);


switch narg
    case 1
        outparam = arg{1};
    case 0
        outparam = 'rhomulambda';
    otherwise
        warning(['wrong number of arguments: ', num2str(narg)]);
end

end