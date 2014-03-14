
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
ux = ux + vx*5*dt;                                                              % IS THIS CORRECT / ACCURATE ENOUGH? --> test with source freq. and dt.
uy = uy + vy*5*dt;                               
uz = uz + vz*5*dt;
% This MAY be wrong: if the displacement at the end of the forward
% calculation hasn't 'died out' yet, my implementation in the code
% artificially sets displacement at T_end at zero (i.e. in
% 'initialise_dynamic_fields.m'. 
% However: I think that the adjoint method requires u and v to be zero at
% T_end anyway, which would mean it's perfectly all right.                      % final displacement possibly needs to be saved

% strain tensor
[duxdx,duxdz,duzdx,duzdz] = grad_v_PSV(ux,uz,dx,dz,nx,nz,order);
% could do this for uy with the function gradient() too.                                
duydx = dx_v(uy,dx,dz,nx,nz,order);
duydz = dz_v(uy,dx,dz,nx,nz,order);


%% get necessary forward fields

% get forward velocity field (remember: saved backwards in time)

vy_fw_snapshot = squeeze(vy_forward(n/5,:,:));
vx_fw_snapshot = squeeze(vx_forward(n/5,:,:));
vz_fw_snapshot = squeeze(vz_forward(n/5,:,:));


% calculate forward displacement field
% PRETTY DAMN INEFFICIENT!!!

ux_fw_snapshot = sum(vx_forward(n/5:end,:,:),1)*5*dt;
uy_fw_snapshot = sum(vy_forward(n/5:end,:,:),1)*5*dt;                           % is it even correct? Do we get the right displacements at all?! Here we're
uz_fw_snapshot = sum(vz_forward(n/5:end,:,:),1)*5*dt;                           % calculating the displacement field from only 1/5 of the velocity field..!


% strain tensor

[duxdx_fw,duxdz_fw,duzdx_fw,duzdz_fw] = grad_v_PSV(ux_fw_snapshot, ...
                                    uz_fw_snapshot,dx,dz,nx,nz,order);
% could do this for uy with the function gradient() too.
duydx_fw = dx_v(uy_fw_snapshot,dx,dz,nx,nz,order);
duydz_fw = dz_v(uy_fw_snapshot,dx,dz,nx,nz,order);                                                                             % HOW does it work for SH?!



%% rho -- density               (calculated from velocities)

% interaction
interaction_vy=vy.*vy_fw_snapshot;
interaction_vx=vx.*vx_fw_snapshot;
interaction_vz=vz.*vz_fw_snapshot;

% kernels
Krho_y=Krho_y-interaction_vy*dt; % The minus sign is needed because we run backwards in time.
Krho_x=Krho_x-interaction_vx*dt;
Krho_z=Krho_z-interaction_vz*dt;

Krho_SH =Krho_y;
Krho_PSV=Krho_x+Krho_z;

Krho = [Krho_x Krho_y Krho_z Krho_PSV Krho_SH];

%% mu -- shear modulus      	(calculated from gradient of displacement)

% interaction
interaction_muPSV = 2*duxdx*duxdx_fw + 2*duzdz*duzdz_fw ...
                    + (duxdz + duzdx) * (duzdx_fw + duxdz_fw);
interaction_muSH  = 1/2 * (duydx*duydx_fw + duydz*duydz_fw);

% kernels
Kmu_PSV = Kmu_PSV + interaction_muPSV*dt;
Kmu_SH = Kmu_SH + interaction_muSH*dt;

Kmu = [Kmu_PSV Kmu_SH];



%% lambda -- lamé's parameter	(calculated from divergence of displacement)
 
% interaction
interaction_lambdaPSV = (duxdx + duzdz) * (duxdx_fw + duzdz_fw);

% kernel
Klambda_PSV = Klambda_PSV + interaction_lambdaPSV*dt;

Klambda = Klambda_PSV;


