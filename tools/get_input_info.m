function [project_name, axrot, apply_hc, use_grav, ...
    use_matfile_startingmodel, starting_model, ...
    parametrisation, param_plot, ...
    rec_g, X, Z, misfit_type, normalise_misfits, stepInit] = get_input_info

% function that gives the project name from the file input_parameters.m

path(path,'../input');
input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

end