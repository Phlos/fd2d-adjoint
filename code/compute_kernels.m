
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
%
% Trying to make it a function (27 Oct 2014), where Kernel should be the
% total kernels only, and K the PSV / SH sub-kernels.
% function [Kernel, K] compute_kernels(vx, vy, vz, v_forward, ...)
%==========================================================================

sfe = store_fw_every;

%% get necessary dynamic fields

%- adjoint ----------------------------------------------------------------

% strain tensor
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
[duxdx,duxdz,duzdx,duzdz] = grad_v_PSV(ux,uz,dx,dz,nx,nz,order);
end
% could do this for uy with the function gradient() too.                                
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
duydx = dx_v(uy,dx,dz,nx,nz,order);
duydz = dz_v(uy,dx,dz,nx,nz,order);
end


%- forward ----------------------------------------------------------------

% forward velocity field (remember: saved backwards in time)
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    vx_fw = squeeze(v_fw.x((n+(sfe-1))/sfe,:,:));
    vz_fw = squeeze(v_fw.z((n+(sfe-1))/sfe,:,:));
%     vx_fw = squeeze(v_fw.x(n/sfe,:,:));
%     vz_fw = squeeze(v_fw.z(n/sfe,:,:));
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    vy_fw = squeeze(v_fw.y((n+(sfe-1))/sfe,:,:));
%     vy_fw = squeeze(v_fw.y(n/sfe,:,:));
end

% forward displacement field
if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    ux_fw = squeeze(u_fw.x((n+(sfe-1))/sfe,:,:));
    uz_fw = squeeze(u_fw.z((n+(sfe-1))/sfe,:,:));
%     ux_fw = squeeze(u_fw.x(n/sfe,:,:));
%     uz_fw = squeeze(u_fw.z(n/sfe,:,:));
end
if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    uy_fw = squeeze(u_fw.y((n+(sfe-1))/sfe,:,:));
%     uy_fw = squeeze(u_fw.y(n/sfe,:,:));
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

% KERNEL CALCULATION:
% generally: K = K - 0.5 * interaction * sfe * dt
% - the factor 0.5 comes from the misfit functional :  1/2 * int ( u - u0 )^2
% - sfe * dt because we are integrating over time and this is the step length
% - the minus sign remains a well-kept secret which I keep confusing myself
%   about each time I look at it.
%      --- Nienke Blom, 23-7-2015

% NOTE about the minus sign:
% all the kernels are calculated with a minus sign, while Andreas' book
% says only rho should be minus. I think this is correct - because now
% the kernels are what I expect them to be. I am still (11-4-2014) 
% utterly confused about this whole minus sign business, but I think 
% the following is the case:
% - the adjoint field is calculated backwards in time, which means that
%   displacements are calculated from velocities in the wrong way.
% - to counteract this, the kernels which make use of the displacement
%   fields (mu and lambda) have to be corrected with an extra minus.
%      --- Nienke Blom, 11-4-2014

    %% rho -- density               (calculated from velocities)
    
    % interaction
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        interaction.rho.y = vy.*vy_fw;
%         interaction.rho.y = vy.*vy_fw * dx * dz;
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        interaction.rho.x = vx.*vx_fw;
        interaction.rho.z = vz.*vz_fw;
%         interaction.rho.x = vx.*vx_fw * dx * dz;
%         interaction.rho.z = vz.*vz_fw * dx * dz;
        interaction.rho.PSV = interaction.rho.x + interaction.rho.z;
    end
    
    % kernels
    
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        K.rho.SH = K.rho.SH - interaction.rho.y*sfe*dt;   
%         K.rho.SH = K.rho.SH - 0.5 * interaction.rho.y*sfe*dt;   
    end
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        K.rho.x   = K.rho.x - interaction.rho.x*sfe*dt;
        K.rho.z   = K.rho.z - interaction.rho.z*sfe*dt;
%         K.rho.x   = K.rho.x - 0.5 * interaction.rho.x*sfe*dt;
%         K.rho.z   = K.rho.z - 0.5 * interaction.rho.z*sfe*dt;
        K.rho.PSV = K.rho.x+K.rho.z;
    end
    
    % previously I had all with +, i.e. the minus of the Krho was minused:
    % Here, two minus signs canceled each other out. Normally, the density
    % kernel would be -v*v_fw (see Andreas' book for details).
    % Normally, v = (u(t+dt) - u(t)) / dt, but we've calculated the adjoint 
    % v backwards in time, so that the executed calculation instead becomes 
    % v_adj = (u(t) - u(t+dt)) / dt. Therefore we need another minus to 
    % compensate.
    % => doesn't this mean that the mu and lambda kernels should have an
    %    extra minus?!! because those are the quantities derived from the v
    %    field which is calculated directly in the wave propagation. Think
    %    about this.
    %      --- Nienke Blom, ~20-3-2014
    
    
    
    %% mu -- shear modulus      	(calculated from gradient of displacement)
    
    % interaction
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        interaction.mu.PSV = ( 2*duxdx.*duxdx_fw + 2*duzdz.*duzdz_fw ...
                            + (duxdz + duzdx) .* (duzdx_fw + duxdz_fw) ); % ...
%                             * dx * dz;
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        interaction.mu.SH  = (1 * (duydx.*duydx_fw + duydz.*duydz_fw) ); % ...
%                             * dx * dz;
    end
    
    % kernels
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        K.mu.PSV = K.mu.PSV - interaction.mu.PSV*sfe*dt;
%         K.mu.PSV = K.mu.PSV - 0.5 * interaction.mu.PSV*sfe*dt;
    end
    if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
        K.mu.SH = K.mu.SH - interaction.mu.SH*sfe*dt;
%         K.mu.SH = K.mu.SH - 0.5 * interaction.mu.SH*sfe*dt;
    end
    
    
    
    
    %% lambda -- lam�'s parameter	(calculated from divergence of displacement)
    
    if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
        % interaction
        interaction.lambda.PSV = ((duxdx + duzdz) .* (duxdx_fw + duzdz_fw)); % ...
%                                 * dx * dz;
        
        % kernel
        K.lambda.PSV = K.lambda.PSV - interaction.lambda.PSV*sfe*dt;
%         K.lambda.PSV = K.lambda.PSV - 0.5 * interaction.lambda.PSV*sfe*dt;
    end
    
%% fill out the kernels which are not calculated but which one may want to plot

if (strcmp(wave_propagation_type,'SH'))
    K.mu.PSV = zeros(size(K.mu.SH));
    K.rho.PSV = zeros(size(K.rho.SH));
    K.lambda.PSV = zeros(size(K.rho.SH));
end

if (strcmp(wave_propagation_type,'PSV'))
    K.mu.SH = zeros(size(K.mu.PSV));
    K.rho.SH = zeros(size(K.rho.PSV));
end

%- total kernels, used for the inversion based on rho mu lambda
%  parametrisation.
K.rho.total = K.rho.PSV + K.rho.SH;
K.mu.total = K.mu.PSV + K.mu.SH;
K.lambda.total = K.lambda.PSV;

