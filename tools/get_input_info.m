function [project_name, axrot, apply_hc, parametrisation, rec_g, X, Z] = get_input_info

% function that gives the project name from the file input_parameters.m


path(path,'../input');
input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

end