
function view_models(modeltype)

% function to view what the models look like

path(path,'../code')

input_parameters;
set_figure_properties;

[mu,rho,lambda]=define_material_parameters(nx,nz,modeltype);               

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

plot_model;

end