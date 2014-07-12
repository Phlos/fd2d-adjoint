function model = calculate_model_properties(rho, axrot)

% function that calculates some overall properties of the model, such as
% mass, moment of inertia, ...
% -- Nienke Blom, 9-7-2014
%
% INPUT:
% - rho - nx * nz double with the model density
%
% OUTPUT:
% - model - struct w/:  model.V volume
%                       model.W r^2 weighted volume (wrt axis of rot.)
%                       model.Z r^4 weighted volume (wrt axis of rot.)
%                       model.M total mass
%                       model.I total moment of inertia (wrt axis of rot.)
%
% NOTE: I could make this prettier with just one if-loop by setting
%       max(r(:)) as the equivalent of Lx or Lz, but that only works if the
%       axis of rotation is indeed x=0 or z=0. If there happen to be
%       addable constants for another axis of rotation for each of these 
%       parameters (such as for the moment of inertia, sth like 'parallel 
%       axis theorem'????), then this would work, but it requires some
%       extra calculations.
%


%- prepare necessary information
path(path,'../input');
path(path,'../tools');
path(path,'../code');
path(path,'../code/propagation');

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

if (axrot == 'x')
    r = Z';
elseif(axrot == 'z')
    r = X';
else
    error('the rotation axis was not recognised')
end

%- calculate volume V
model.V = Lx * Lz;


%- calculate   W = int-over-volume[ r^2 ]d3x (r^2 weighted volume)
model.W = sum(sum(r .* r));
%        and   Z = int-over-volume[ r^4 ]d3x (r^4 weighted volume)
model.Z = sum(sum(r .* r .* r .* r));

% if (axrot == 'x')
%     model.W = 1/3 * Lx * Lz^3;
%     model.Z = 1/5 * Lx * Lz^5;
% elseif(axrot == 'z')
%     model.W = 1/3 * Lz * Lx^3;
%     model.Z = 1/5 * Lz * Lx^5;
% else
%     error('the rotation axis was not recognised')
% end


%- calculate model total mass M
model.M = sum(sum(rho(2:end,2:end)))*dx*dz;

%- calculate model total moment of inertia I

dI = rho .* r .* r;
% sizedI = size(dI)
% sizerho= size(rho)
model.I = sum(sum(rho .* r .* r));

end