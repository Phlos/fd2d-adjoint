
function fig_mod = view_models(modeltype)

% function to view what the predefined models look like.
%
% INPUT:
% - modeltype:	number defining the model
%
% OUTPUT:
% - a figure with the model parameters plotted
% - fig_mod:    the figure handle of the figure containing the model
% 

path(path,'../code')

input_parameters;
set_figure_properties_bothmachines;

[model.mu,model.rho,model.lambda]=define_material_parameters(nx,nz,modeltype);               

% [X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

fig_mod = plot_model(model);

end