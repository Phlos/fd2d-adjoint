function [fig_mod, Model_diff] = plot_model_diff(Model1, Model2, varargin)

% function that plots the difference model between Model 1 and Model 2
% Model_diff = Model1 - Model2

if isnumeric(Model1)
    Model1 = update_model(Model1);
end
if isnumeric(Model2)
    Model2 = update_model(Model2);
end

% varargin{:}
[outparam, plot_UM_separate] = checkargs(varargin);

switch outparam
    case 'rhomulambda'
        Model_diff.rho = Model1.rho - Model2.rho;
        Model_diff.mu = Model1.mu - Model2.mu;
        Model_diff.lambda = Model1.lambda - Model2.lambda;
    case 'rhovsvp'
        Model1 = change_parametrisation('rhomulambda','rhovsvp',Model1);
        Model2 = change_parametrisation('rhomulambda','rhovsvp',Model2);
        Model_diff.rho = Model1.rho - Model2.rho;
        Model_diff.vs = Model1.vs - Model2.vs;
        Model_diff.vp = Model1.vp - Model2.vp;
    otherwise
        error('parametrisation not recognised');
end


% plotting the model with middle = 0 meaning that the white of the plot is
% the value zero.
middle.rho    = 0;
middle.mu     = 0;
middle.lambda = 0;
fig_mod = plot_model(Model_diff,middle,varargin{:});


end




function [outparam, plot_UM_separate] = checkargs(arg)

% defines the output parametrisation based on whether any is given (default
% is 'rhomulambda')

narg = length(arg);

% default values
outparam = 'rhomulambda';
plot_UM_separate = false;

for ii = 1:narg
    if ischar(arg{ii})
        outparam = arg{ii};
    elseif islogical(arg{ii})
        plot_UM_separate = arg{ii};
    else
        error('input argument not recognised')
    end
end


% switch narg
%     case 1
%         outparam = arg{1};
%     case 0
%         outparam = 'rhomulambda';
%     otherwise
%         warning(['wrong number of arguments: ', num2str(narg)]);
% end

end