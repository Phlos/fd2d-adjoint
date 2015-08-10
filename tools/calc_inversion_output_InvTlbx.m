function [InvProps] = calc_inversion_output_InvTlbx(iter, currentIterFields, prevIterFields, InvProps, K_total, Model, jm, gm, varargin)

% calculate some inversion properties to be plotted in order to monitor the
% inversion development.


input_parameters; 
Model_real = checkargs(varargin(:));
Ktotal_iter = K_total(iter);
Kseis_iter = currentIterFields.Kseis;
Kg_iter = currentIterFields.Kg;
if iter > 1
    Ktotal_prev = K_total(iter-1);
    Kseis_prev = prevIterFields.Kseis;
    Kg_prev     = prevIterFields.Kg;
end

%% MISFITS

InvProps.misfit(iter) = jm;
InvProps.misfitseis(iter) = currentIterFields.misfit_seis;
InvProps.misfitgrav(iter) = currentIterFields.misfit_grav;

InvProps.gradmag(iter) = norm(gm);


%% MODEL NORMS
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


%% KERNEL NORM
% L2 norm of current kernels:
% if iter < InvProps.niter
    InvProps.norm.Kseis(iter) = norm(currentIterFields.Kseis.rho.total(:)) + ...
                               norm(currentIterFields.Kseis.mu.total(:)) + ...
                               norm(currentIterFields.Kseis.lambda.total(:));
    InvProps.norm.Kg(iter) = norm(currentIterFields.Kg(:));
%     InvProps.norm.Kseis(iter)  = NaN;
%     InvProps.norm.Kg(iter)     = NaN;
    if strcmp(parametrisation, 'rhomulambda');
        InvProps.norm.Ktotal(iter) = norm(Ktotal_iter.rho.total(:)) + ...
            norm(Ktotal_iter.mu.total(:)) + ...
            norm(Ktotal_iter.lambda.total(:));
    elseif strcmp(parametrisation, 'rhovsvp');
        InvProps.norm.Ktotal(iter) = norm(Ktotal_iter.rho2.total(:)) + ...
            norm(Ktotal_iter.vs2.total(:)) + ...
            norm(Ktotal_iter.vp2.total(:));
    else
        error('calc_inversion_output_InvTlbx: parametrisation of knls not recognised');
    end
% else
%     InvProps.norm.Kseis(iter) = NaN;
%     InvProps.norm.Kg(iter) = NaN;
%     InvProps.norm.Ktotal(iter) = NaN;
% end
    

%% KERNEL ANGLE
    %- calculate angles between model updates:

% angle between kernels
if (length(K_total) > 1 && iter > 1)
    par = fieldnames(K_total);
    KnlTot_now = [Ktotal_iter.(par{1}).total Ktotal_iter.(par{2}).total Ktotal_iter.(par{3}).total];
    KnlTot_prv = [Ktotal_prev.(par{1}).total Ktotal_prev.(par{2}).total Ktotal_prev.(par{3}).total];
    KnlS_now = [Kseis_iter.(par{1}).total Kseis_iter.(par{2}).total Kseis_iter.(par{3}).total];
    KnlS_prv = [Kseis_prev.(par{1}).total Kseis_prev.(par{2}).total Kseis_prev.(par{3}).total];
    
    totNtimesP                  = KnlTot_prv(:)' * KnlTot_now(:);
    totNorm                     = norm(KnlTot_prv(:))*norm(KnlTot_now(:));
    InvProps.angle.Ktotal(iter) = acos(totNtimesP / totNorm);
    
    if strcmp(use_seis,'yesseis')
        seisNtimesP                 = KnlS_prv(:)' * KnlS_now(:);
        seisNorm                    = norm(KnlS_prv(:))*norm(KnlS_now(:));
        InvProps.angle.Kseis(iter)  = acos(seisNtimesP / seisNorm);
    else
        InvProps.angle.Kseis(iter)  = NaN;
    end
    
    if strcmp(use_grav,'yes');
        gravNtimesP                 = Kg_prev(:)' * Kg_iter(:);
        gravNorm                    = norm(Kg_prev(:))*norm(Kg_iter(:));
        InvProps.angle.Kg(iter)     = acos(gravNtimesP / gravNorm);
    else
        InvProps.angle.Kg(iter)     = NaN;
    end
else
    InvProps.angle.Ktotal(iter) = NaN;
    InvProps.angle.Kseis(iter)  = NaN;
    InvProps.angle.Kg(iter)  = NaN;
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