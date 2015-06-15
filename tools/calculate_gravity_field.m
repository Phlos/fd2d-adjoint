

function [g, fig_grav] = calculate_gravity_field(rho, rec_grav, varargin)

% calculate gravity field
%
% INPUT
%
% - receiver information: rec_grav.x,.z for all receivers.
% - body information: density (as a function of x, z) [height width dx dz 
%   all derived from input_parameters
%
% OUTPUT
% - gravity vector at each receiver
% 
% TODO:
% - change rec_grav to optional input.


% plot gravity field or not? Default: yes.
plotornot = plot_or_not(varargin(:));

%- prepare necessary information
% path(path,'./input');
% path(path,'./tools');
% path(path,'./code');
% path(path,'./code/propagation');

% gravitational constant
% 6.67384 * 10-11 m3 kg-1 s-2
G = 6.67384e-11;

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

nrec = size(rec_grav.x,2);

% calculate dm: mass for each block of model
dm = rho * dx * dz;

for i = 1:nrec
%     disp(['receiver ', num2str(i), ' out of ', num2str(nrec)]);
    % calculate distance vector r{i}.x, r{i}.z ??
    r{i}.x = rec_grav.x(i) - X';
    r{i}.z = rec_grav.z(i) - Z';
    
    % calculate length or the vectors r for each point
    r{i}.length = sqrt(r{i}.x.^2 + r{i}.z.^2);
    
    % THIS IS FOR GRAVITY DUE TO A 'SHEET' IN THE X-Z PLANE!!
    % calculate dg{i}.x,z: gravity increment for each block of model for each receiver
    dg{i}.x = - dm ./ r{i}.length .^3 .* r{i}.x * G;
    dg{i}.z = - dm ./ r{i}.length .^3 .* r{i}.z * G;
    
%     % FOR GRAVITY WITH Y STRETCHING TO INFINITY AT BOTH SIDES
%     % calculate dg{i}.x,z: gravity increment for each block of model for each receiver
%     dg{i}.x = - dm ./ r{i}.length .^3 .* r{i}.x * G    .* 2 .* r{i}.length;
%     dg{i}.z = - dm ./ r{i}.length .^3 .* r{i}.z * G    .* 2 .* r{i}.length;

    % calculate g(i): gravity for each receiver
%     g{i}.x = sum(dg{i}.x(:));
%     g{i}.z = sum(dg{i}.z(:));
% disp 'calculating total gravity field'
    g.x(i) = sum(dg{i}.x(:));
    g.z(i) = sum(dg{i}.z(:));
    g.mag(i) = sqrt(g.x(i)^2 + g.z(i)^2);
    
end

% disp(['maximum gravity accelleration: ', num2str(max(g.mag)), ' m/s^2'])

% create dummy gravity field to be able to plot the gravity field 
% (the plot function now has static input but really I should change the
% function plot_gravity_quivers so that you can either plot one field,
% or the difference between two (or even multiple fields, maybe?)
dum.x = zeros(size(g.x));
dum.z = zeros(size(g.z));
dum.mag = zeros(size(g.mag));

if strcmp(plotornot, 'yesplot')
    fig_grav = plot_gravity_quivers(rec_grav, g, dum, X, Z, rho);
else
    fig_grav = NaN;
end

end

function plotornot = plot_or_not(args)

% plot gravity field or not? Default: yes.
plottext = {'yesplot','noplot'};

if isempty(args)
    plotornot = 'yesplot';
elseif (length(args) == 1 &&  ischar(args{1}) && any(strcmp(args{1}, plottext)) )
    plotornot = args{1};
else
    error('plotornot variable input not recognised');
end

end