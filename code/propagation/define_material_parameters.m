function [mu,rho,lambda]=define_material_parameters(nx,nz,model_type)

%==========================================================================
% generate material parameters lambda, mu [N/m^2] and rho [kg/m^3]
%
% input:    grid points of the velocity and density field in x-direction (nx) and z-direction (nz)
%           model_type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%==========================================================================

if (model_type==1)
    
    % homogeneous
    rho=3000.0*ones(nx,nz);
    mu=4.8e10*ones(nx,nz);
    lambda=mu;
    
elseif (model_type==2)
    
    % homogeneous with density perturbation
    rho=3000.0*ones(nx,nz);
    mu=4.8e10*ones(nx,nz);
    lambda=mu;
    
    rho(98:102,123:127)=rho(98:102,123:127)+2000.0;
    
elseif (model_type==3)
    
    % two halves with different rho, mu, lambda
    rho=3000.0*ones(nx,nz);
    mu=2.8e10*ones(nx,nz);
    
    rho(1:round(nx/2),:)=rho(1:round(nx/2),:)+200.0;
    mu(1:round(nx/2),:)=mu(1:round(nx/2),:)+3.5e10;
    lambda=mu;
    
elseif (model_type==4)
    
    rho=3000.0*ones(nx,nz);
    mu=2.8e10*ones(nx,nz);
    
    rho(1:round(nx/2),:)=rho(1:round(nx/2),:)+200.0;
    mu(1:round(nx/2),:)=mu(1:round(nx/2),:)+3.5e10;
    lambda=mu;
    
    rho(98:102,123:127)=rho(98:102,123:127)+2000.0;
    
elseif (model_type==5)
    
    rho=3000.0*ones(nx,nz);
    mu=2.0e10*ones(nx,nz);
    lambda=mu;
    
    for k=100:150
        mu(:,k)=mu(:,k)+(k-100)*10.0e8;
    end
    for k=151:nz
        mu(:,k)=mu(:,150);
    end
   
elseif (model_type==6)
    
    rho=3000.0*ones(nx,nz);
    mu=2.0e10*ones(nx,nz);
    lambda=mu;
    
    for k=100:150
        mu(:,k)=mu(:,k)+(k-100)*10.0e8;
    end
    for k=151:nz
        mu(:,k)=mu(:,150);
    end
    
    rho(98:102,123:127)=rho(98:102,123:127)+2000.0;
    
elseif (model_type==11)
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
elseif (model_type==12)
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    left = round(nx/2-nx/20);
    right = round(nx/2+nx/20);
    top = round(nz/2+nz/20);
    bottom = round(nz/2-nz/20);
    
    mu(left:right,bottom:top) = mu(left:right,bottom:top) + 1e10;
    
elseif (model_type==13)
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    left = round(nx/2-nx/20);
    right = round(nx/2+nx/20);
    top = round(nz/2+nz/20);
    bottom = round(nz/2-nz/20);
    
    rho(left:right,bottom:top) = rho(left:right,bottom:top) + 1e3;
    
elseif (model_type==100)
    
    rho=3000.0*ones(nx,nz);
    mu=ones(nx,nz);
%     lambda=3e10*ones(nx,nz);
    
    mu(1:330,:)=3.675e10;
    mu(331:end,:)=2.7e10;
    lambda=mu;
    
else
    
    load(['models/mu_' model_type]);
    load(['models/rho_' model_type]);
    
end


