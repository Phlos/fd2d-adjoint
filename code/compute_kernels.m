
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
% Klambda:          lamé parameter kernel:  _PSV
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



%% get necessary adjoint fields from velocity field

% displacement field (calc'd every 5th timestep)
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
ux = ux + vx*5*dt;                                                              % IS THIS CORRECT / ACCURATE ENOUGH? --> test with source freq. and dt.
uz = uz + vz*5*dt;
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
uy = uy + vy*5*dt;  
end

% This MAY be wrong: if the displacement at the end of the forward
% calculation hasn't 'died out' yet, my implementation in the code
% artificially sets displacement at T_end at zero (i.e. in
% 'initialise_dynamic_fields.m'. 
% However: I think that the adjoint method requires u and v to be zero at
% T_end anyway, which would mean it's perfectly all right.                      % final displacement possibly needs to be saved

% strain tensor
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
[duxdx,duxdz,duzdx,duzdz] = grad_v_PSV(ux,uz,dx,dz,nx,nz,order);
end
% could do this for uy with the function gradient() too.                                
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
duydx = dx_v(uy,dx,dz,nx,nz,order);
duydz = dz_v(uy,dx,dz,nx,nz,order);
end

%% get necessary forward fields

% get forward velocity field (remember: saved backwards in time)
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
vy_fw_snapshot = squeeze(vy_forward(n/5,:,:));
end
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
vx_fw_snapshot = squeeze(vx_forward(n/5,:,:));
vz_fw_snapshot = squeeze(vz_forward(n/5,:,:));
end

% calculate forward displacement field
% PRETTY DAMN INEFFICIENT!!!

if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
ux_fw_snapshot = squeeze(sum(vx_forward(n/5:end,:,:),1)*5*dt);
uz_fw_snapshot = squeeze(sum(vz_forward(n/5:end,:,:),1)*5*dt);                           % calculating the displacement field from only 1/5 of the velocity field..!
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
uy_fw_snapshot = squeeze(sum(vy_forward(n/5:end,:,:),1)*5*dt);                           % is it even correct? Do we get the right displacements at all?! Here we're
end

% strain tensor

if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
[duxdx_fw,duxdz_fw,duzdx_fw,duzdz_fw] = grad_v_PSV(ux_fw_snapshot, ...
                                    uz_fw_snapshot,dx,dz,nx,nz,order);
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
% could do this for uy with the function gradient() too.
duydx_fw = dx_v(uy_fw_snapshot,dx,dz,nx,nz,order);
duydz_fw = dz_v(uy_fw_snapshot,dx,dz,nx,nz,order);                                                                             % HOW does it work for SH?!
end


%% calculate travel time kernels
if strcmp(kerneltype,'traveltime')
    %% rho -- density               (calculated from velocities)
    
    % interaction
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    interaction_vy=vy.*vy_fw_snapshot;
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    interaction_vx=vx.*vx_fw_snapshot;
    interaction_vz=vz.*vz_fw_snapshot;
    end
    
    % kernels
    
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    K.rho.SH=K.rho.SH+interaction_vy*5*dt; % The minus sign is needed because we run backwards in time.
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    K.rho.x=K.rho.x+interaction_vx*5*dt;
    K.rho.z=K.rho.z+interaction_vz*5*dt;
    % K.rho.SH =K.rho.y;
    K.rho.PSV=K.rho.x+K.rho.z;
    end
    

    
    
    
    %% mu -- shear modulus      	(calculated from gradient of displacement)
    
    % interaction
    
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    interaction_muPSV = 2*duxdx.*duxdx_fw + 2*duzdz.*duzdz_fw ...
        + (duxdz + duzdx) .* (duzdx_fw + duxdz_fw);
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    interaction_muSH  = 1/2 * (duydx.*duydx_fw + duydz.*duydz_fw);
    end
    
    % kernels
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    K.mu.PSV = K.mu.PSV + interaction_muPSV*5*dt;
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    K.mu.SH = K.mu.SH + interaction_muSH*5*dt;
    end
    % K.mu.y = K.mu.SH;
    
    
    
    
    %% lambda -- lamé's parameter	(calculated from divergence of displacement)
    
    % interaction
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    interaction_lambdaPSV = (duxdx + duzdz) .* (duxdx_fw + duzdz_fw);
    
    % kernel
    K.lambda.PSV = K.lambda.PSV + interaction_lambdaPSV*5*dt;
    end
else
    error('Sorry, you have to specify you want travel time kernels')
    
end
