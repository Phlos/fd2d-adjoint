function [InvProps] = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model)


        InvProps.misfitseis(iter) = InvProps.misfit_seis{iter}.normd;
        InvProps.misfitgrav(iter) = InvProps.misfit_g(iter).normd;

    
    %- calculate angles between model updates:
%     for ie = 2:length(Kg)
    if (length(Kg) > 1 && iter < InvProps.niter)
        kernelskeerelkaar = Kg{iter-1}(:)'*Kg{iter}(:);
        normz = norm(Kg{iter-1}(:))*norm(Kg{iter}(:));
        InvProps.angle.Kg(iter) = acos(kernelskeerelkaar / normz);
    else
        InvProps.angle.Kg(iter) = NaN;
    end
    
%     for iter = 1:length(Kseis)
    if (length(Kseis) > 1 && iter < InvProps.niter)
        
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
    
    if (length(K_total) > 1 && iter < InvProps.niter)
        Ktransf{iter} = change_parametrisation_kernels('rhomulambda','rhovsvp',K_total(iter),Model(iter));
        totaalkernel{iter} = [Ktransf{iter}.rho2.total Ktransf{iter}.vs2.total Ktransf{iter}.vp2.total];
%     end
%     for iter = 2:length(K_total)
        kernelskeerelkaar = totaalkernel{iter-1}(:)'*totaalkernel{iter}(:);
        normz = norm(totaalkernel{iter-1}(:))*norm(totaalkernel{iter}(:));
        InvProps.angle.Ktotal(iter) = acos(kernelskeerelkaar / normz);
    else
        InvProps.angle.Ktotal(iter) = NaN;
    end
    
end