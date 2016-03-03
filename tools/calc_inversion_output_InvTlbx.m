function [InvProps] = calc_inversion_output_InvTlbx(iter, currentIterFields, prevIterFields, InvProps, K_total, Model, jm, gm, varargin)

% calculate some inversion properties to be plotted in order to monitor the
% inversion development.


input_parameters;
Model_real = checkargs(varargin(:));
Ktotal_iter = K_total(iter);
if strcmp(use_seis, 'yesseis')
    Kseis_iter = currentIterFields.Kseis;
else
    Kseis_iter = NaN;
end
if strcmp(use_grav, 'yes')
    Kg_iter = currentIterFields.Kg;
else
    Kg_iter = NaN;
end
if iter > 1
    Ktotal_prev = K_total(iter-1);
    if strcmp(use_seis, 'yesseis')
        Kseis_prev = prevIterFields.Kseis;
    else
        Kseis_prev = NaN;
    end
    if strcmp(use_grav, 'yes')
        Kg_prev     = prevIterFields.Kg;
    else
        Kg_prev = NaN;
    end
    
end

if strcmp(parametrisation, 'rhovsvp');
    Ktotal_iter = change_parametrisation_kernels('rhomulambda', 'rhovsvp', Ktotal_iter, Model_bg);
    Ktotal_prev = change_parametrisation_kernels('rhomulambda', 'rhovsvp', Ktotal_prev, Model_bg);
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



if(isstruct(Model_real))
    
    % models / abbrev:
    MR = Model_real;
    Mi = Model(iter);
    MB = update_model(bg_model_type); %Model_bg;
    MRR = change_parametrisation('rhomulambda', 'rhovsvp', MR);
    MiR = change_parametrisation('rhomulambda', 'rhovsvp', Mi );
    MBR = change_parametrisation('rhomulambda', 'rhovsvp', MB);
    
    % L2 norm (abs): [model(i) - model_real] / model_real
    % param rho-mu-lambda
    InvProps.modeldifnFromTruePar.rho(iter) = norm( (Mi.rho(:) - MR.rho(:)) ./ MR.rho(:) );
    InvProps.modeldifnFromTruePar.mu(iter) = norm( (Mi.mu(:)  - MR.mu(:))  ./ MR.mu(:) );
    InvProps.modeldifnFromTruePar.lambda(iter) = norm( (Mi.lambda(:) - MR.lambda(:)) ./ MR.lambda(:) );
    % param rho-vs-vp
    InvProps.modeldifnFromTruePar.vs(iter) = norm( (MiR.vs(:)  - MRR.vs(:))  ./ MRR.vs(:) );
    InvProps.modeldifnFromTruePar.vp(iter) = norm( (MiR.vp(:) - MRR.vp(:)) ./ MRR.vp(:) );
    % total diff
    InvProps.modeldifnFromTrue(iter) = InvProps.modeldifnFromTruePar.rho(iter) + ...
        InvProps.modeldifnFromTruePar.mu(iter) + InvProps.modeldifnFromTruePar.lambda(iter);
    
    % L2 norm (new): |[model(i) - model_real]| / |[model_real - model_start]| )
    L2nnew.rho    = norm( (Mi.rho(:) - MR.rho(:)) , 2)  ./ norm( (MR.rho(:) - MB.rho(:))  , 2);
    L2nnew.mu     = norm( (Mi.mu(:) - MR.mu(:)) , 2)  ./ norm( (MR.mu(:) - MB.mu(:))  , 2);
    L2nnew.lambda = norm( (Mi.lambda(:) - MR.lambda(:)) , 2)  ./ norm( (MR.lambda(:) - MB.lambda(:))  , 2);
    L2nnew.vs     = norm( (MiR.vs(:) - MRR.vs(:)) , 2)  ./ norm( (MRR.vs(:) - MBR.vs(:))  , 2);
    L2nnew.vp     = norm( (MiR.vp(:) - MRR.vp(:)) , 2)  ./ norm( (MRR.vp(:) - MBR.vp(:))  , 2);
    L2nnew.total  = L2nnew.rho + L2nnew.vs + L2nnew.vp;
else
    
    % L2 norm (old)
    % param rho-mu-lambda
    InvProps.modeldifnFromTruePar.rho(iter) = NaN;
    InvProps.modeldifnFromTruePar.mu(iter) = NaN;
    InvProps.modeldifnFromTruePar.lambda(iter) = NaN;
    % param rho-vs-vp
    InvProps.modeldifnFromTruePar.vs(iter) = NaN;
    InvProps.modeldifnFromTruePar.vp(iter) = NaN;
    % total
    InvProps.modeldifnFromTrue(iter) = NaN;
    
     % L2 norm new
    L2nnew.rho    = NaN;
    L2nnew.mu     = NaN;
    L2nnew.lambda = NaN;
    L2nnew.vs     = NaN;
    L2nnew.vs     = NaN;
end

InvProps.L2norm_normd.rho(iter)    = L2nnew.rho;
InvProps.L2norm_normd.mu(iter)     = L2nnew.mu;
InvProps.L2norm_normd.lambda(iter) = L2nnew.lambda;
InvProps.L2norm_normd.vs(iter)     = L2nnew.vs;
InvProps.L2norm_normd.vp(iter)     = L2nnew.vp;
InvProps.L2norm_normd.total_rvv(iter)  = (L2nnew.rho + L2nnew.vs + L2nnew.vp)/3;
InvProps.L2norm_normd.total_rml(iter)  = (L2nnew.rho + L2nnew.mu + L2nnew.lambda)/3;



%% KERNEL NORM
% L2 norm of current kernels:
% if iter < InvProps.niter
if strcmp(use_seis, 'yesseis');
    InvProps.norm.Kseis(iter) = norm(Kseis_iter.rho.total(:)) + ...
                               norm(Kseis_iter.mu.total(:)) + ...
                               norm(Kseis_iter.lambda.total(:));
else
    InvProps.norm.Kseis(iter) = NaN;
end
if strcmp(use_grav, 'yes')
    InvProps.norm.Kg(iter) = norm(Kg_iter);
else
    InvProps.norm.Kg(iter) = NaN;
end
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
    

%% KERNEL ANGLE
    %- calculate angles between model updates:

% angle between kernels
if (length(K_total) > 1 && iter > 1)
    par = fieldnames(K_total);
    KnlTot_now = [Ktotal_iter.(par{1}).total Ktotal_iter.(par{2}).total Ktotal_iter.(par{3}).total];
    KnlTot_prv = [Ktotal_prev.(par{1}).total Ktotal_prev.(par{2}).total Ktotal_prev.(par{3}).total];
    if strcmp(use_seis,'yesseis')
        KnlS_now = [Kseis_iter.(par{1}).total Kseis_iter.(par{2}).total Kseis_iter.(par{3}).total];
        KnlS_prv = [Kseis_prev.(par{1}).total Kseis_prev.(par{2}).total Kseis_prev.(par{3}).total];
    end
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