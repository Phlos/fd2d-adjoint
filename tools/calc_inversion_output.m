function [InvProps] = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model, varargin)

Model_real = checkargs(varargin(:));

% calculate some inversion properties to be plotted in order to monitor the
% inversion development.

input_parameters; 
if iter < InvProps.niter;
    Kg_iter = filter_kernels(Kg{iter},smoothgwid);
    Ktotal_iter.rho.total = filter_kernels(K_total(iter).rho.total,smoothgwid);
    Ktotal_iter.mu.total = filter_kernels(K_total(iter).mu.total,smoothgwid);
    Ktotal_iter.lambda.total = filter_kernels(K_total(iter).lambda.total,smoothgwid);
end


% L2 norm( [model(i) - model(i-1)] / model(1) )
if iter>1
    InvProps.modeldifn(iter) =   norm( (Model(iter).rho(:) - Model(iter-1).rho(:)) ./ Model(1).rho(:) ) ...
        + norm( (Model(iter).mu(:)  - Model(iter-1).mu(:))  ./ Model(1).mu(:) ) ...
        + norm( (Model(iter).lambda(:) - Model(iter-1).lambda(:)) ./ Model(1).lambda(:) );

else
    InvProps.modeldifn(iter) = NaN;
end

% L2 norm( [model(i) - model_real] / model_real )
if(isstruct(Model_real))
    InvProps.modeldifnFromTrue(iter) = norm( (Model(iter).rho(:) - Model_real.rho(:)) ./ Model_real.rho(:) ) ...
        + norm( (Model(iter).mu(:)  - Model_real.mu(:))  ./ Model_real.mu(:) ) ...
        + norm( (Model(iter).lambda(:) - Model_real.lambda(:)) ./ Model_real.lambda(:) );
else
    InvProps.modeldifnFromTrue(iter) = NaN;
end

% L2 norm of current kernels:
if iter < InvProps.niter
InvProps.norm.Kseis(iter) = norm(Kseis(iter).rho.total(:)) + ...
                           norm(Kseis(iter).mu.total(:)) + ...
                           norm(Kseis(iter).lambda.total(:));
InvProps.norm.Kg(iter) = norm(Kg_iter(:));
InvProps.norm.Ktotal(iter) = norm(Ktotal_iter.rho.total(:)) + ...
                            norm(Ktotal_iter.mu.total(:)) + ...
                            norm(Ktotal_iter.lambda.total(:));
else
    InvProps.norm.Kseis(iter) = NaN;
    InvProps.norm.Kg(iter) = NaN;
    InvProps.norm.Ktotal(iter) = NaN;
end
    
    %- calculate angles between model updates:
%     angle between gravity kernels
if strcmp(use_grav,'yes')
    if (length(Kg) > 1 && iter > 1 && iter < InvProps.niter)
        kernelskeerelkaar = Kg{iter-1}(:)'*Kg{iter}(:);
        normz = norm(Kg{iter-1}(:))*norm(Kg{iter}(:));
        InvProps.angle.Kg(iter) = acos(kernelskeerelkaar / normz);
    else
        InvProps.angle.Kg(iter) = NaN;
    end
end
    
%    angle between seismic kernels
if (length(Kseis) > 1 && iter > 1 && iter < InvProps.niter)
    
    Ktransf{iter} = change_parametrisation_kernels('rhomulambda','rhovsvp',Kseis(iter),Model(iter));
    Ktransf{iter-1} = change_parametrisation_kernels('rhomulambda','rhovsvp',Kseis(iter-1),Model(iter-1));
    totaalkernel{iter} = [Ktransf{iter}.rho2.total Ktransf{iter}.vs2.total Ktransf{iter}.vp2.total];
    totaalkernel{iter-1} = [Ktransf{iter-1}.rho2.total Ktransf{iter-1}.vs2.total Ktransf{iter-1}.vp2.total];
    
    kernelskeerelkaar = totaalkernel{iter-1}(:)'*totaalkernel{iter}(:);
    normz = norm(totaalkernel{iter-1}(:))*norm(totaalkernel{iter}(:));
    InvProps.angle.Kseis(iter) = acos(kernelskeerelkaar / normz);
else
    InvProps.angle.Kseis(iter) = NaN;
end

% angle between total kernels
if (length(K_total) > 1 && iter > 1 && iter < InvProps.niter)
    Ktransf{iter} = change_parametrisation_kernels('rhomulambda','rhovsvp',K_total(iter),Model(iter));
    totaalkernel{iter} = [Ktransf{iter}.rho2.total Ktransf{iter}.vs2.total Ktransf{iter}.vp2.total];
    
    kernelskeerelkaar = totaalkernel{iter-1}(:)'*totaalkernel{iter}(:);
    normz = norm(totaalkernel{iter-1}(:))*norm(totaalkernel{iter}(:));
    InvProps.angle.Ktotal(iter) = acos(kernelskeerelkaar / normz);
else
    InvProps.angle.Ktotal(iter) = NaN;
end
    
end

function Mod_real = checkargs(args)

narg = numel(args);

if narg == 0
    Mod_real = NaN;
elseif narg == 1 && isstruct(args{1})
    Mod_real = args{1};
end

end