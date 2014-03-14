% Initialises all the dynamic fields (in terms of zeros(nx,nz) and such)

% Do I really need to make these ifs? Does it hurt to initialise dynamic
% fiels which may not be used at all? -- NAB 13-3-2014

% velocity fields (for one timestep)
if(strcmp(wave_propagation_type,'SH'))
    vy=zeros(nx,nz);
    % displacement
    uy=zeros(nx,nz);
elseif(strcmp(wave_propagation_type,'PSV'))
    vx=zeros(nx,nz);
    vz=zeros(nx,nz);
    % displacement
    ux=zeros(nx,nz);
    uz=zeros(nx,nz);
elseif(strcmp(wave_propagation_type,'both'))
    vy=zeros(nx,nz);
    vx=zeros(nx,nz);
    vz=zeros(nx,nz);
    % displacement
    ux=zeros(nx,nz);
    uy=zeros(nx,nz);
    uz=zeros(nx,nz);
end

if(strcmp(simulation_mode,'forward'))
    % stored velocity fields (a nt/5 by nx by nz matrix)
    if(strcmp(wave_propagation_type,'SH'))
        vy_forward=zeros(nt/5,nx,nz);
    elseif(strcmp(wave_propagation_type,'PSV'))
        vx_forward=zeros(nt/5,nx,nz);
        vz_forward=zeros(nt/5,nx,nz);
    elseif(strcmp(wave_propagation_type,'both'))
        vy_forward=zeros(nt/5,nx,nz);
        vx_forward=zeros(nt/5,nx,nz);
        vz_forward=zeros(nt/5,nx,nz);
    end
end

% stress fields
if(strcmp(wave_propagation_type,'SH'))
    sxy=zeros(nx-1,nz); % stress component sigma_xy (nx-1 because of stag grid)
    szy=zeros(nx,nz-1); % stress component sigma_zy
elseif(strcmp(wave_propagation_type,'PSV'))
    sxx=zeros(nx,nz);                                                            %%%%%%%%% LET OP hier heb ik ze aangepast
    szz=zeros(nx,nz);                                                            % sxx: nx-1 --> nx;  szz nz-1 --> nz
    sxz=zeros(nx,nz);         
elseif(strcmp(wave_propagation_type,'both'))
    sxy=zeros(nx-1,nz); % stress component sigma_xy (nx-1 because of stag grid)
    szy=zeros(nx,nz-1); % stress component sigma_zy
    sxx=zeros(nx,nz);                                                            %%%%%%%%% LET OP hier heb ik ze aangepast
    szz=zeros(nx,nz);                                                            % sxx: nx-1 --> nx;  szz nz-1 --> nz
    sxz=zeros(nx,nz);  
end