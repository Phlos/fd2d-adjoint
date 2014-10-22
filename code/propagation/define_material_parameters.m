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
    
elseif (model_type==10) % homogeneous like Tromp et al 2005
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
elseif (model_type==11) % tromp05-homogeneous + rect. mu anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    % rectangular mu anomaly
    left = round(nx/2-nx/20);
    right = round(nx/2+nx/20);
    top = round(nz/2+nz/20);
    bottom = round(nz/2-nz/20);
    
    mu(left:right,bottom:top) = mu(left:right,bottom:top) + 1e10;
    
elseif (model_type==12) % tromp05-homogeneous + rect. rho anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    % rectangular rho anomaly
    left = round(nx/2-nx/20);
    right = round(nx/2+nx/20);
    top = round(nz/2+nz/20);
    bottom = round(nz/2-nz/20);
    
    rho(left:right,bottom:top) = rho(left:right,bottom:top) + 1e3;
    
elseif (model_type==14) % gaussian central rho anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    gwid = round(0.05 * max([nx nz]));
    filt = fspecial('gaussian',[nx nz],gwid);
    filt2 = filt / max(filt(:));
    rho = rho + filt2 * 1.0e3;
    
elseif (model_type==15) % gaussian central mu anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    gwid = round(0.05 * max([nx nz]));
    filt = fspecial('gaussian',[nx nz],gwid);
    filt2 = filt / max(filt(:));
    mu = mu + filt2 * 1.0e10;
    
elseif (model_type==17) % gaussian off-central rho anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    % location of the centre of the anomaly
    anom.dxperc=0.40; % anomaly x distance from the lower left corner
                  % as a fraction of the X length
    anom.dzperc=0.42; % same for z distance
    % height & width of the anomaly
    anom.dx = round(anom.dxperc*nx);
    anom.dz  = round(anom.dzperc*nz);
    
    % make anomaly
    gwid = round(0.05 * max([nx nz])); % standard deviation of gaussian:
                                       % size of the gaussian blur of the 
                                       % rho anomaly, fraction of the total
                                       % domain size
    filt.width = 8*gwid; % 8*gwid is enough to taper out visibly
    filt.height = filt.width; % don't see why it should ever be different..
    filt.f = fspecial('gaussian',[filt.width filt.height],gwid); 
    filt.fnorm = filt.f / max(filt.f(:)); % normalised filter (max value is 1)
    
    % add the anomaly to the density field at (anom.dx, anom.dz)
    rho = add_crop_matrix(rho, filt.fnorm*10^3, anom.dx, anom.dz);
    
    
elseif (model_type==18) % gaussian off-central mu anomaly
    
    % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    dxz=0.15;
    
    gwid = round(0.05 * max([nx nz]));
    filt = fspecial('gaussian',[floor((1-dxz)*nx) floor((1-dxz)*nz)],gwid);
    filt2 = filt / max(filt(:));
%     size_filt2 = size(filt2)
%     size_rhocut = size(rho(ceil((1-dxz)*nx)+1:end, ceil((1-dxz)*nz)+1:end))
    mu(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) = ...
    mu(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) + filt2 * 1.0e10;

elseif (model_type==21) % gaussian off-central rho_v anomaly
%% VP-VS-RHO parametrisation

    % Tromp et al, 2005, rho-mu-lambda
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
        % rho-vs-vp
        vp    = sqrt((lambda + 2*mu) ./ rho);
        vs    = sqrt(mu ./ rho);
        rho2 = rho;
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
%     % adding the gaussian rho_v anomaly
%     dxz=0.15;
%     gwid = round(0.05 * max([nx nz]));
%     filt = fspecial('gaussian',[floor((1-dxz)*nx) floor((1-dxz)*nz)],gwid);
%     filt2 = filt / max(filt(:));
%     rho2(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) = ...
%     rho2(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) + filt2 * 1.0e3;

    % location of the centre of the anomaly
    anom.dxperc=0.40; % anomaly x distance from the lower left corner
                  % as a fraction of the X length
    anom.dzperc=0.42; % same for z distance
    % height & width of the anomaly
    anom.dx = round(anom.dxperc*nx);
    anom.dz  = round(anom.dzperc*nz);

    % make anomaly
    gwid = round(0.05 * max([nx nz])); % standard deviation of gaussian:
                                       % size of the gaussian blur of the 
                                       % rho anomaly, fraction of the total
                                       % domain size
    filt.width = 8*gwid; % 8*gwid is enough to taper out visibly
    filt.height = filt.width; % don't see why it should ever be different..
    filt.f = fspecial('gaussian',[filt.width filt.height],gwid); 
    filt.fnorm = filt.f / max(filt.f(:)); % normalised filter (max value is 1)
    
    % add the anomaly to the density field at (anom.dx, anom.dz)
    rho2 = add_crop_matrix(rho2, filt.fnorm*10^3, anom.dx, anom.dz);

    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==100) % layered: left = high velocity, right = low vel.
    
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


