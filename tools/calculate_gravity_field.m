

function [r, r_len, dm, dg, g] = calculate_gravity_field(rec_grav, rho)

% calculate gravity field
%
% INPUT
%
% - body information: height width dx dz density (as a function of x, z,
%   indices) --> comes from input_parameters
% - receiver information: rec_g_x, rec_g_z for all receivers.
%
% OUTPUT
% - gravity vector at each receiver
% 

%- prepare necessary information
path(path,'../input');
path(path,'../tools');
path(path,'../code');
path(path,'../code/propagation');

% gravitational constant
% 6.67384 × 10-11 m3 kg-1 s-2
G = 6.67384e-11;

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

nrec = size(rec_grav.x,2);

% calculate dm: mass for each block of model
dm = rho * dx * dz;

for i = 1:nrec
    disp(['receiver ', num2str(i), ' out of ', num2str(nrec)]);
    % calculate distance vector r{i}.x, r{i}.z ??
    r{i}.x = rec_grav.x(i) - X';
    r{i}.z = rec_grav.z(i) - Z';
    
    
    % calculate length or the vectors r for each point
    r_len{i} = sqrt(r{i}.x.^2 + r{i}.z.^2);
    
    % calculate dg{i}(iks,zet,comp): gravity increment for each block of model for each receiver
    disp 'calculating dg x'
    dg{i}.x = - dm ./ r_len{i}.^3 .* r{i}.x * G;
    disp 'calculating dg z'
    dg{i}.z = - dm ./ r_len{i}.^3 .* r{i}.z * G;

    % calculate g(i): gravity for each receiver
%     g{i}.x = sum(dg{i}.x(:));
%     g{i}.z = sum(dg{i}.z(:));
disp 'calculating total gravity field'
    g.x(i) = sum(dg{i}.x(:));
    g.z(i) = sum(dg{i}.z(:));
    g.mag(i) = sqrt(g.x(i)^2 + g.z(i)^2);
    
end




end