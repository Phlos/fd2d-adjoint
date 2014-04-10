
% calculate sensitivity kernels or misfit kernels for all manner of things.
%
% 'input'           (i.e. quantities that need to be present)
% -------
% vx, vy, vz:       adjoint velocity field in current timestep
% v{xyz}_forward:   forward velocity field (from stored matfile)
% dx, dz, dt:       parameters governing the gradient & derivative calcs
%
% output:
% -------
% Krho:             density kernels:        _x _y _z _PSV _SH
% Kmu:              shear modulus kernels:  _PSV _SH
% Klambda:          lamÃ© parameter kernel:  _PSV
%
% Note the following:
% Note that kernels for density can be calculated for x, y and z dirs
% separately. For lambda and mu, however, those of x and z seismograms are
% coupled so you need both x and z measurements to get the full kernel. For
% lambda, there is no sensitivity in shear waves at all (remember: vs =
% sqrt( mu/rho ), so no sensitivity in the y component of seismograms!
%
% -- Nienke Blom, 13 March 2014
%==========================================================================



%% get necessary dynamic fields



% adjoint

% strain tensor
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
[duxdx,duxdz,duzdx,duzdz] = grad_v_PSV(ux,uz,dx,dz,nx,nz,order);
end
% could do this for uy with the function gradient() too.                                
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
duydx = dx_v(uy,dx,dz,nx,nz,order);
duydz = dz_v(uy,dx,dz,nx,nz,order);
end


% forward

% forward velocity field (remember: saved backwards in time)
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    vx_fw = squeeze(v_fw.x(n/5,:,:));
    vz_fw = squeeze(v_fw.z(n/5,:,:));
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    vy_fw = squeeze(v_fw.y(n/5,:,:));
end

% forward displacement field
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    ux_fw = squeeze(u_fw.x(n/5,:,:));
    uz_fw = squeeze(u_fw.z(n/5,:,:));
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    uy_fw = squeeze(u_fw.y(n/5,:,:));
end

% forward strain tensor
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    [duxdx_fw,duxdz_fw,duzdx_fw,duzdz_fw] = grad_v_PSV(ux_fw, ...
                                                  uz_fw,dx,dz,nx,nz,order);
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    duydx_fw = dx_v(uy_fw,dx,dz,nx,nz,order);   % could do this for uy with the function gradient() too.
    duydz_fw = dz_v(uy_fw,dx,dz,nx,nz,order);                                                                             % HOW does it work for SH?!
end






%% calculate kernels

% if strcmp(kerneltype,'traveltime')
    
    %% rho -- density               (calculated from velocities)
    
    % interaction
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        interaction.rho.y = vy.*vy_fw;
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        interaction.rho.x = vx.*vx_fw;
        interaction.rho.z = vz.*vz_fw;
        interaction.rho.PSV = interaction.rho.x + interaction.rho.z;
    end
    
    % kernels
    
    % NOTE:
    % Here, two minus signs cancel each other out. Normally, the density
    % kernel would be -v*v_fw (see Andreas' book for details).
    % Normally, v = (u(t+dt) - u(t)) / dt, but we've calculated the adjoint 
    % v backwards in time, so that the executed calculation instead becomes 
    % v = (u(t) - u(t+dt)) / dt. Therefore we need another minus to 
    % compensate.
    % => doesn't this mean that the mu and lambda kernels should have an
    %    extra minus?!! because those are the quantities derived from the v
    %    field which is calculated directly in the wave propagation. Think
    %    about this.
    
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        K.rho.SH = K.rho.SH + interaction.rho.y*5*dt;   
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        K.rho.x   = K.rho.x + interaction.rho.x*5*dt;
        K.rho.z   = K.rho.z + interaction.rho.z*5*dt;
        K.rho.PSV = K.rho.x+K.rho.z;
    end
    

    
    
    
    %% mu -- shear modulus      	(calculated from gradient of displacement)
    
    % interaction
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        interaction.mu.PSV = 2*duxdx.*duxdx_fw + 2*duzdz.*duzdz_fw ...
                            + (duxdz + duzdx) .* (duzdx_fw + duxdz_fw);
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        interaction.mu.SH  = 1 * (duydx.*duydx_fw + duydz.*duydz_fw);
    end
    
    % kernels
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        K.mu.PSV = K.mu.PSV + interaction.mu.PSV*5*dt;
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        K.mu.SH = K.mu.SH + interaction.mu.SH*5*dt;
    end
    
    
    
    
    %% lambda -- lamé's parameter	(calculated from divergence of displacement)
    
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        % interaction
        interaction.lambda.PSV = (duxdx + duzdz) .* (duxdx_fw + duzdz_fw);
        
        % kernel
        K.lambda.PSV = K.lambda.PSV + interaction.lambda.PSV*5*dt;
    end
    
    
    
    
    
% else
%     error('Sorry, you have to specify you want travel time kernels')
%     
% end
