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
    gwid = round(0.05 * max([nx nz])); % size of the anomaly
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
    
    % make anomaly
    gwid = round(0.05 * max([nx nz])); % size of the anomaly
    filt.width = 8*gwid; % 8*gwid is enough to taper out visibly
    filt.height = filt.width; % don't see why it should ever be different..
    filt.f = fspecial('gaussian',[filt.width filt.height],gwid); 
    filt.fnorm = filt.f / max(filt.f(:)); % normalised filter (max value is 1)
    
    % add the anomaly to the density field at (anom.dx, anom.dz)
    mu = add_crop_matrix(mu, filt.fnorm*10^3, anom.dx, anom.dz);
    
%     gwid = round(0.05 * max([nx nz]));
%     filt = fspecial('gaussian',[floor((1-dxz)*nx) floor((1-dxz)*nz)],gwid);
%     filt2 = filt / max(filt(:));
% %     size_filt2 = size(filt2)
% %     size_rhocut = size(rho(ceil((1-dxz)*nx)+1:end, ceil((1-dxz)*nz)+1:end))
%     mu(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) = ...
%     mu(ceil(dxz*nx)+1:end, ceil(dxz*nz)+1:end) + filt2 * 1.0e10;

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
    
    % location of the centre of the anomaly
    anom.dxperc=0.40; % anomaly x distance from the lower left corner
                  % as a fraction of the X length
    anom.dzperc=0.42; % same for z distance
    % height & width of the anomaly
    anom.dx = round(anom.dxperc*nx);
    anom.dz  = round(anom.dzperc*nz);

    % make anomaly
    gwid = round(0.05 * max([nx nz])); % size of the anomaly
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
    
elseif (model_type==31) % five 'rand' rho2 anomalies (rho2 = rho in rho-vs-vp)

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
    
%     willekeurig.x = [ 0.7577    0.7431    0.3922    0.6555    0.1712];
    willekeurig.x = [ 0.1576    0.9706    0.4854    0.9572    0.8003];
    willekeurig.z = [ 0.8235    0.6948    0.3171    0.9502    0.0344];
    
    for i = 1:size(willekeurig.x,2)
%         disp(['i = ', num2str(i), ' out of ']);
        anom{i}.dxperc = willekeurig.x(i);
        anom{i}.dzperc = willekeurig.z(i);
        
        
        % height & width of the anomaly
        anom{i}.dx = round(anom{i}.dxperc*nx);
        anom{i}.dz  = round(anom{i}.dzperc*nz);

        % make anomaly
        gwid = round(0.05 * max([nx nz])); % size of the anomaly
        filt.width = 8*gwid; % 8*gwid is enough to taper out visibly
        filt.height = filt.width; % don't see why it should ever be different..
        filt.f = fspecial('gaussian',[filt.width filt.height],gwid);
        filt.fnorm = filt.f / max(filt.f(:)); % normalised filter (max value is 1)

        % add the anomaly to the density field at (anom.dx, anom.dz)
        rho2 = add_crop_matrix(rho2, filt.fnorm*10^3, anom{i}.dx, anom{i}.dz);

    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==41) % ten 'rand' rho2 anomalies (rho2 = rho in rho-vs-vp)
                        % 5x positive and 5x negative

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
    
    rho2 = add_10randanoms(rho2, 1e3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==50) % PREM background model
                        % model values will be sampled at height above CMB!
                        % IMPORTANT: 
                        % so don't make the model higher than 2891 km!!
                        
    [rho_new, vs_new, vp_new] = load_PREM();    
    
    %- convert to rho,mu,lambda
    rho     = rho_new;
    mu      = vs_new .^ 2 .* rho_new;
    lambda  = rho_new .* ( vp_new.^2 - 2* vs_new.^2);
    
    
    
elseif (model_type==51) % PREM background model + 10 rand +&- rho2 anoms
                        % model values will be sampled at height above CMB!
                        % IMPORTANT: 
                        % so don't make the model higher than 2891 km!!
                      
    [rho2, vs, vp] = load_PREM();
    rho2 = add_10randanoms(rho2, 1e3);
   
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==52) % PREM background model + 10 rand +&- vs anoms
                        % model values will be sampled at height above CMB!
                        % IMPORTANT: 
                        % so don't make the model higher than 2891 km!!
                        
    [rho2, vs, vp] = load_PREM();
    vs = add_10randanoms(vs, 1e3);
  
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==53) % PREM background model + 10 SMALL rand +&- vs anoms
                        % model values will be sampled at height above CMB!
                        % IMPORTANT: 
                        % so don't make the model higher than 2891 km!!

    [rho2, vs, vp] = load_PREM();
    vs = add_10randanoms(vs, 0.01*max(vs(:)));
    
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
    
elseif (model_type == 101) % test model with tiny rho anomaly
    
        % Tromp et al, 2005
    rho    = 2600*ones(nx,nz);     % kg/m3
    mu     = 2.66e10*ones(nx,nz);  % Pa
    lambda = 3.42e10*ones(nx,nz);  % Pa
    % => vp = 5797.87759 m/s
    % => vs = 3198.55736 m/s
    
    % rho-vs-vp
    vp    = sqrt((lambda + 2*mu) ./ rho);
    vs    = sqrt(mu ./ rho);
    rho2 = rho;
    
    rho2(98:102,123:127)=rho2(98:102,123:127)+50.0;
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    

elseif (model_type == 102) % Evangelos: ring shaped model
    
    vp_ring = 5000;
    vs_ring  = 3000; % m/s
    rho_ring = 2600; % kg/m3
    
    vp_out  = 5000;
    vs_out  = 1; % m/s
    rho_out = 2600; % kg/m3
    
    % rho-vs-vp
    vp    = vp_out*ones(nx,nz);
%     vs    = sqrt(mu ./ rho);
    vs   = vs_out*ones(nx,nz);
%     rho2 = zeros(nx,nz);
    rho2 = rho_out*ones(nx,nz);
    
    % outer circle
    c_outer.mid.x = nx/2;
    c_outer.mid.z = nz/2;
    c_outer.radius = min(nx,nz)/2;
    
    % inner circle
    c_inner.mid.x = nx/2;
    c_inner.mid.z = nz/2;
    c_inner.radius = min(nx/2,nz/2)/2;
    
    for ii = 1:size(vs,1)
        for jj = 1:size(vs,2)
            location = sqrt((ii-c_outer.mid.x)^2 + (jj-c_outer.mid.z)^2);
            if location < c_outer.radius
                vp(ii,jj) = vp_ring;
                vs(ii,jj) = vs_ring;
                rho2(ii,jj) = rho_ring;
                if location < c_inner.radius
                    vp(ii,jj) = vp_out;
                    vs(ii,jj) = vs_out;
                    rho2(ii,jj) = rho_out;
                end
            end
        end
    end
    
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
else
    
    load(['models/mu_' model_type]);
    load(['models/rho_' model_type]);
    
end

end

function [rho2, vs, vp] = load_PREM()


% NOTE: PREM read from table with columns like
    % http://ds.iris.edu/ds/products/emc-prem/PREM_1s.csv
    % units on that website are in g/cm^3, m/s, m/s
    
    %- make a grid using the normal Lx, Lz, nx, nz
    input_parameters;
    [X,Z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
    
    %- load PREM data from file
    PREM = csvread('~/Dropbox/DensityInversion/PREM-reference-model/PREM_1s_nowater.csv');
    
    %- depth coordinate manipulation
    depth = PREM(:,2);      % in km
    % change double occurrences of depth so that interpolation doesn't flip
    dif_dep = find([NaN; diff(depth)] <= 0);
    depth(dif_dep) = depth(dif_dep)+1E-5;
    %- convert depth to height above CMB
    h_CMB = 2891 - depth;   % in km!

    %- make a grid using PREM values (nonuniform grid)
    %  (1000* for units km -> m , g/cm^3 -> kg/m^3 , km/s -> m/s)
    [X_PREM, Z_PREM] = meshgrid( X(1,:),1000* h_CMB(:) );
    [~, rho_PREM] = meshgrid(X(1,:),1000* PREM(:,3));
    [~, vp_PREM]  = meshgrid(X(1,:),1000* PREM(:,4));     % in reality vp_vertical 
    [~, vs_PREM]  = meshgrid(X(1,:),1000* PREM(:,6));     % in reality vs_vertical
    
        
    %- interpolate PREM vp, vs, rho to our grid sampling vp, vs, rho
    rho_new = interp2(X_PREM, Z_PREM, rho_PREM, X,Z);
    vp_new  = interp2(X_PREM, Z_PREM, vp_PREM,  X,Z);
    vs_new  = interp2(X_PREM, Z_PREM, vs_PREM,  X,Z);
    
    % transpose
    rho2 = rho_new'; vp = vp_new'; vs = vs_new';
    
end

function mod_out = add_10randanoms(mod, anommax)

nx = size(mod,1);
nz = size(mod,2);

% random locations of anomalies
    willekeurig.x = [0.906812, 0.493549, 0.499576, 0.175434, 0.441924, ...
                     0.668359, 0.541729, 0.444079, 0.348643, 0.927618];
    willekeurig.z = [0.261980, 0.930472, 0.645794, 0.534028, 0.907064, ...
                     0.631545, 0.039304, 0.194939, 0.657794, 0.749388];
    
    % add random anomalies to model
    for ii = 1:size(willekeurig.x,2)
        % anomaly strength
%         anom{i}.strength = 0.10 * 2600;
        anom{ii}.strength = anommax;
        
        anom{ii}.dxperc = willekeurig.x(ii);
        anom{ii}.dzperc = willekeurig.z(ii);
       
        % x and z location of the anomaly
        anom{ii}.dx = round(anom{ii}.dxperc*nx);
        anom{ii}.dz  = round(anom{ii}.dzperc*nz);

        % make anomaly
        gwid = round(0.05 * max([nx nz])); % size of the anomaly
        filt.width = 8*gwid; % 8*gwid is enough to taper out visibly
        filt.height = filt.width; % don't see why it should ever be different..
        filt.f = fspecial('gaussian',[filt.width filt.height],gwid);
        filt.fnorm = anom{ii}.strength * filt.f / max(filt.f(:)); % normalised filter (max value is 1)
        % make half of the anomalies negative
        if (rem(ii,2)==0)
            filt.fnorm = - filt.fnorm;
        else
        end
        
        % add the anomaly to the VS field at (anom.dx, anom.dz)
        mod = add_crop_matrix(mod, filt.fnorm, anom{ii}.dx, anom{ii}.dz);
    end
    
    mod_out = mod;

end