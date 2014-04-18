% run the actual wavefield propagation


%% initialise dynamic fields ----------------------------------------------
% both forward and adjoint
initialise_dynamic_fields;  % this just makes all dynamic field (v, stress,
                            % derivatives of v and stress wherever needed
                            % with zeros(dimensions).


%% initialise absorbing boundary taper a la Cerjan ------------------------

absbound=ones(nx,nz);
init_absbound;

%%
%==========================================================================
% iterate
%==========================================================================

disp 'iterating...'

position_figures;


%%
for n=1:nt
    
    %- compute divergence of current stress tensor ------------------------
    %  forward -- unit of force: [kg m^-2 s^-2] = [N] 
    
    if(strcmp(wave_propagation_type,'SH'))
        DSY=div_s(sxy,szy,dx,dz,nx,nz,order);
    elseif(strcmp(wave_propagation_type,'PSV'))
        [DSX,DSZ]=div_s_PSV(sxx,szz,sxz,dx,dz,nx,nz,order);
    elseif(strcmp(wave_propagation_type,'both'))
        DSY=div_s(sxy,szy,dx,dz,nx,nz,order);
        [DSX,DSZ]=div_s_PSV(sxx,szz,sxz,dx,dz,nx,nz,order);
    end
    
%       test.DSX_before_stf1(n) = DSX(src_x_id(1),src_z_id(1));
%       test.DSX_before_stf2(n) = DSX(src_x_id(2),src_z_id(2));
    
    %- add point sources --------------------------------------------------
    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green') ||  strcmp(simulation_mode,'adjoint') )
    
        for i=1:ns
            
            if(strcmp(wave_propagation_type,'SH'))
                DSY(src_x_id(i),src_z_id(i))= DSY(src_x_id(i),src_z_id(i)) + stf(2,i,n);
            elseif(strcmp(wave_propagation_type,'PSV'))
                DSX(src_x_id(i),src_z_id(i))= DSX(src_x_id(i),src_z_id(i)) + stf(1,i,n);                   %%%%%%%%% Beetje krukkig zo... kan het mooier?
                DSZ(src_x_id(i),src_z_id(i))= DSZ(src_x_id(i),src_z_id(i)) + stf(3,i,n);
            elseif(strcmp(wave_propagation_type,'both'))
                DSY(src_x_id(i),src_z_id(i))= DSY(src_x_id(i),src_z_id(i)) + stf(2,i,n);
                DSX(src_x_id(i),src_z_id(i))= DSX(src_x_id(i),src_z_id(i)) + stf(1,i,n);                   %%%%%%%%% Beetje krukkig zo... kan het mooier?
                DSZ(src_x_id(i),src_z_id(i))= DSZ(src_x_id(i),src_z_id(i)) + stf(3,i,n);
            end

        end
        

    end
    
    
    %- update velocity field ----------------------------------------------
    
%     test.vx_source1 = vx(src_x_id(1),src_z_id(1)) + ...
%             dt*DSX(src_x_id(1),src_z_id(1)) / rho(src_x_id(1),src_z_id(1));
%     test.vx_source2 = vx(src_x_id(2),src_z_id(2)) + ...
%             dt*DSX(src_x_id(2),src_z_id(2)) / rho(src_x_id(2),src_z_id(2));
% size(test_vx_source)
    
    if(strcmp(wave_propagation_type,'SH'))
        vy=vy+dt*DSY./rho;
    elseif(strcmp(wave_propagation_type,'PSV'))
        vx=vx+dt*DSX./rho;
        vz=vz+dt*DSZ./rho;
    elseif(strcmp(wave_propagation_type,'both'))
        vy=vy+dt*DSY./rho;
        vx=vx+dt*DSX./rho;
        vz=vz+dt*DSZ./rho;
    end
    
                        % just a test to see whether DS and v produce values
                        % [max(max(DS)),max(max(DSX)),max(max(DSZ)); max(max(v)),max(max(vx)),max(max(vz))]
    
%     test.stf(n) = stf(1,1,n);
%     test.DSX1(n) = DSX(src_x_id(1),src_z_id(1));
%     test.DSX2(n) = DSX(src_x_id(2),src_z_id(2));
%     test.vx(n)  = vx(src_x_id(1),src_z_id(1));
%     disp(['stf ',num2str(stf(1,1,n)), ...
%               '  DSX  ', num2str(DSX(src_x_id(1),src_z_id(1))) ...
%               '  vx   ', num2str(vx(src_x_id(1),src_z_id(1)))]) % ...
%               '  testvx ', num2str(test_vx_source)])
    
    %- apply absorbing boundary taper -------------------------------------
    
    if(strcmp(wave_propagation_type,'SH'))
        vy=vy.*absbound;
    elseif(strcmp(wave_propagation_type,'PSV'))
        vx=vx.*absbound;
        vz=vz.*absbound;
    elseif(strcmp(wave_propagation_type,'both'))
        vy=vy.*absbound;
        vx=vx.*absbound;
        vz=vz.*absbound;
    end
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    if(strcmp(wave_propagation_type,'SH'))
        sxy=sxy+ dt* mu.* dx_v(vy,dx,dz,nx,nz,order);
        szy=szy+ dt* mu.* dz_v(vy,dx,dz,nx,nz,order);
    elseif(strcmp(wave_propagation_type,'PSV'))
        % calc gradient of velocity field
        [dvxdx,dvxdz,dvzdx,dvzdz]=grad_v_PSV(vx,vz,dx,dz,nx,nz,order);
        % stress update
        sxx=sxx+dt*( (lambda+2*mu).*dvxdx(:,:) + lambda.*dvzdz(:,:) );
        szz=szz+dt*( (lambda+2*mu).*dvzdz(:,:) + lambda.*dvxdx(:,:) );
        sxz=sxz+dt*( mu.*(dvxdz(:,:) + dvzdx(:,:)) );
    elseif(strcmp(wave_propagation_type,'both'))
        sxy=sxy + dt*mu.* dx_v(vy,dx,dz,nx,nz,order);
        szy=szy + dt*mu.* dz_v(vy,dx,dz,nx,nz,order);
        % calc gradient of velocity field (P-SV)
        [dvxdx,dvxdz,dvzdx,dvzdz]=grad_v_PSV(vx,vz,dx,dz,nx,nz,order);
        % stress update (P-SV)
        sxx=sxx+dt*( (lambda+2*mu).*dvxdx(:,:) + lambda.*dvzdz(:,:) );
        szz=szz+dt*( (lambda+2*mu).*dvzdz(:,:) + lambda.*dvxdx(:,:) );
        sxz=sxz+dt*( mu.* (dvxdz(:,:) + dvzdx(:,:)) );
    end
    
    %- compute displacement field -----------------------------------------
    
    if(strcmp(wave_propagation_type,'SH'))
        uy = uy + vy.*dt;
    elseif(strcmp(wave_propagation_type,'PSV'))
        ux = ux + vx.*dt;
        uz = uz + vz.*dt;
    elseif(strcmp(wave_propagation_type,'both'))
        ux = ux + vx.*dt;
        uy = uy + vy.*dt;
        uz = uz + vz.*dt;
    end
    
    
    %% forward calculations
    if (strcmp(simulation_mode,'forward'))
        
        
        % record velocity seismograms -------------------------------------
        
        for k=1:n_receivers
            if(strcmp(wave_propagation_type,'SH'))
                v_rec.y(k,n)=vy(rec_x_id(k),rec_z_id(k));
            elseif(strcmp(wave_propagation_type,'PSV'))
                v_rec.x(k,n)=vx(rec_x_id(k),rec_z_id(k));
                v_rec.z(k,n)=vz(rec_x_id(k),rec_z_id(k));
            elseif(strcmp(wave_propagation_type,'both'))
                v_rec.y(k,n)=vy(rec_x_id(k),rec_z_id(k));
                v_rec.x(k,n)=vx(rec_x_id(k),rec_z_id(k));
                v_rec.z(k,n)=vz(rec_x_id(k),rec_z_id(k));
            end
        end
        
        
        % store time-reversed history ----------------------------------------
        
        % save every 5th velocity&displacement field to the big-ass 3 direction matrix
        if (mod(n,5)==0)
            if(strcmp(wave_propagation_type,'SH'))
                vy_forward(nt/5+1-n/5,:,:)=vy(:,:);
                % displacement
                uy_forward(nt/5+1-n/5,:,:)=uy(:,:);
            elseif(strcmp(wave_propagation_type,'PSV'))
                vx_forward(nt/5+1-n/5,:,:)=vx(:,:);
                vz_forward(nt/5+1-n/5,:,:)=vz(:,:);
                % displacement
                ux_forward(nt/5+1-n/5,:,:)=ux(:,:);
                uz_forward(nt/5+1-n/5,:,:)=uz(:,:);
            elseif(strcmp(wave_propagation_type,'both'))
                vy_forward(nt/5+1-n/5,:,:)=vy(:,:);
                vx_forward(nt/5+1-n/5,:,:)=vx(:,:);
                vz_forward(nt/5+1-n/5,:,:)=vz(:,:);
                % displacement
                uy_forward(nt/5+1-n/5,:,:)=uy(:,:);
                ux_forward(nt/5+1-n/5,:,:)=ux(:,:);
                uz_forward(nt/5+1-n/5,:,:)=uz(:,:);
            end
        end

    
        % plot velocity field every so manyth time step -------------------
        
        if (mod(n,plot_every)==0)
%             disp(['time: ', num2str(n*dt)]);
            plot_velocity_field;
        end
        
        
        
    %% adjoint calculations    
    elseif(strcmp(simulation_mode,'adjoint'))
        
        % plot and compute kernel every 5th iteration -----------------------
        if (mod(n,5)==0)
                        
            % calculate kernels
            compute_kernels;
            
            % plot adjoint fields (i.e. adjoint field + the above)
            if (mod(n,plot_every)==0)
                plot_adjoint_field;
            end
            
        end
    end
    
end