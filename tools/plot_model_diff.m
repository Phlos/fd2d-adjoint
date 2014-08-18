function [Model_diff, fig_mod] = plot_model_diff(Model1, Model2, varargin)

% function that plots the difference model between Model 1 and Model 2
% Model_diff = Model1 - Model2


Model_diff.mu = Model1.mu - Model2.mu;
Model_diff.rho = Model1.rho - Model2.rho;
Model_diff.lambda = Model1.lambda - Model2.lambda;

% max(Model_diff.rho(:))

% plotting the model with the '0' meaning that the white of the plot is
% the value zero.

fig_mod = plot_model(Model_diff,0,varargin{:});

% disp 'To get a correct diff plot, edit the colour scale properties'
% disp 'manually so that the middle value is zero.'

end