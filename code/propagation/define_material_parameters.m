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
    
%     rho(98:102,123:127)=rho(98:102,123:127)+2000.0;
    rho(floor(nx/2 - nx/20):ceil(nx/2 + nx/20) , floor(nz/2 - nz/20):ceil(nz/2 + nz/20) ) = ...
        rho(floor(nx/2 - nx/20):ceil(nx/2 + nx/20) , floor(nz/2 - nz/20):ceil(nz/2 + nz/20) ) + 2000.0;
    
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
    
    anom.dx=0.15;
    anom.dz=0.15;
    
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
    anom.dzperc=0.60; % same for z distance
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
    
    rho2 = add_10randanoms(rho2, 1e3, 1);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
%% PREM bg models    
elseif (model_type==50) % PREM background model
                        % IMPORTANT: 
                        % model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
                        
    [rho_new, vs_new, vp_new] = load_PREM();    
    
    %- convert to rho,mu,lambda
    rho     = rho_new;
    mu      = vs_new .^ 2 .* rho_new;
    lambda  = rho_new .* ( vp_new.^2 - 2* vs_new.^2);
    
    
    
elseif (model_type==51) % PREM background model + 10 rand +&- rho2 anoms
                        % IMPORTANT: 
                        % model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
                      
    [rho2, vs, vp] = load_PREM();
    rho2 = add_10randanoms(rho2, 1000,1);
   
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==52) % PREM background model + 10 rand +&-1000 m/s vs anoms
                        % IMPORTANT: 
                        % model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
                        
    [rho2, vs, vp] = load_PREM();
    vs = add_10randanoms(vs, 1e3,1);
  
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==53) % PREM background model + 10 1% rand +&- vs anoms
                        % IMPORTANT: 
                        % model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!

    [rho2, vs, vp] = load_PREM();
    vs = add_10randanoms(vs, 0.01*max(vs(:)),1);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    
elseif (model_type==54) % PREM background model + 10 1% rand +&- rho2 anoms
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!

    [rho2, vs, vp] = load_PREM();
    rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);

elseif (model_type==55) % PREM background model + 10 1000 kg/m3 rand +&- RHO0 anoms
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!

    [rho2, vs, vp] = load_PREM();
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    rho = add_10randanoms(rho, 1000,1);
    

elseif model_type == 56
    
    [rho2, vs, vp] = load_PREM();
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
    % adding square rho anomaly
     rho(floor(nx/2 - nx/20):ceil(nx/2 + nx/20) , floor(nz/2 - nz/20):ceil(nz/2 + nz/20) ) = ...
        rho(floor(nx/2 - nx/20):ceil(nx/2 + nx/20) , floor(nz/2 - nz/20):ceil(nz/2 + nz/20) ) + 2000.0;
    
    
elseif (model_type==60) % PREM background model + random rho2 AND vs AND vp
                        % (anomaly strength = 1% of largest amount)
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
                        
    [rho2, vs, vp] = load_PREM();
    rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    vs = add_10randanoms(vs, 0.01*max(vs(:)),2);
    vp = add_10randanoms(vp, 0.01*max(vp(:)),3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);

    
    
elseif (model_type==61) % PREM background model + random ONLY vs & vp
                        % (anomaly strength = 1% of largest amount)
                        % NOTE: 
                        % this model can be used as a starting model agains
                        % true model 60, with varying rho2 vs vp.
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
                        
    anom_strength = 0.01;
    [rho2, vs, vp] = load_PREM();
%     rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    vs = add_10randanoms(vs, anom_strength*max(vs(:)),2);
    vp = add_10randanoms(vp, anom_strength*max(vp(:)),3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);

    
elseif (model_type==62) % PREM background model + random ONLY vs & vp
                        % (anomaly strength = 0.5% of largest amount)
                        % NOTE: 
                        % this model can be used as a starting model agains
                        % true model 60, with varying rho2 vs vp.
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    anom_strength = 0.005;
    [rho2, vs, vp] = load_PREM();
%     rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    vs = add_10randanoms(vs, anom_strength*max(vs(:)),2);
    vp = add_10randanoms(vp, anom_strength*max(vp(:)),3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    


elseif (model_type==63) % PREM background model + random ONLY vs & vp
                        % (anomaly strength = 0.5% of largest amount)
                        % NOTE: 
                        % this model can be used as a starting model agains
                        % true model 60, with varying rho2 vs vp.
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    anom_strength = 0.009;
    [rho2, vs, vp] = load_PREM();
%     rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    vs = add_10randanoms(vs, anom_strength*max(vs(:)),2);
    vp = add_10randanoms(vp, anom_strength*max(vp(:)),3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    

elseif (model_type==64) % PREM background model + random ONLY vs & vp
                        % (anomaly strength = 0.5% of largest amount)
                        % NOTE: 
                        % this model can be used as a starting model agains
                        % true model 60, with varying rho2 vs vp.
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    anom_strength = 0.0095;
    [rho2, vs, vp] = load_PREM();
%     rho2 = add_10randanoms(rho2, 0.01*max(rho2(:)),1);
    vs = add_10randanoms(vs, anom_strength*max(vs(:)),2);
    vp = add_10randanoms(vp, anom_strength*max(vp(:)),3);
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
        
    
elseif (model_type==70) % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (6+1), nz / (2+1)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.4] .* spacing; % [to the right, up];
    %- divide into 6+1 horizontally
    locX = spacing(1) * [1:11] + shift(1);
    %- divide into 2/3?+1 vertically
    locZ = spacing(2) * [1:3] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.35 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    rhopos(2).mid = [locX(1), locZ(3)];
    rhoneg(2).mid = [locX(4), locZ(1)];
    rhopos(3).mid = [locX(4), locZ(2)];
    rhoneg(3).mid = [locX(4), locZ(3)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    vspos(2).mid = [locX(2), locZ(3)];
    vsneg(2).mid = [locX(5), locZ(1)];
    vspos(3).mid = [locX(5), locZ(2)];
    vsneg(3).mid = [locX(5), locZ(3)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    vppos(2).mid = [locX(3), locZ(3)];
    vpneg(2).mid = [locX(6), locZ(1)];
    vppos(3).mid = [locX(6), locZ(2)];
    vpneg(3).mid = [locX(6), locZ(3)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2) , ...
                            norm([ii,jj] - rhopos(2).mid , 2)     , ...
                            norm([ii,jj] - rhopos(3).mid , 2)]       );
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2) , ...
                            norm([ii,jj] - rhoneg(2).mid , 2)     , ...
                            norm([ii,jj] - rhoneg(3).mid , 2)]       );
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2) , ...
                            norm([ii,jj] - vspos(2).mid , 2)    , ...
                            norm([ii,jj] - vspos(3).mid , 2) ]      );
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2) , ...
                            norm([ii,jj] - vsneg(2).mid , 2)    , ...
                            norm([ii,jj] - vsneg(3).mid , 2) ]      );
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) , ...
                            norm([ii,jj] - vppos(2).mid , 2)    , ...
                            norm([ii,jj] - vppos(3).mid , 2) ]      );
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2) , ...
                            norm([ii,jj] - vpneg(2).mid , 2)    , ...
                            norm([ii,jj] - vpneg(3).mid , 2) ]      );
            % rho2
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    

    
elseif (model_type==80) % CIRCLES:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+1), nz / (2)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.5] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.35 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2)]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            % rho2
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    
elseif (model_type==81) % CIRCLES only rho:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+1), nz / (2)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.5] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.35 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2)]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            % rho2
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end

        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    
    
elseif (model_type==82) % CIRCLES - like 80 but smaller:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+1), nz / (2)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.5] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.25 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2)]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            % rho2
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    
elseif (model_type==83) % CIRCLES (small type) only rho:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+1), nz / (2)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.5] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.25 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2)]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            % rho2
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end

        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
        
    
elseif (model_type==84) % CIRCLES (small type) only vs vp:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+1), nz / (2)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0 -0.5] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.25 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2)]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);    
    
    
elseif (model_type==85) % CIRCLES and small UM circles:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+2), nz / (2+1)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0.5 -0.2] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.30 * min(nx/(3+1), nz/(2+1));
%     anom.radius = 0.25 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    anomUM.strength = anom.strength;
    anomUM.radius = 0.5 * anom.radius;
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
    rhopos(1).mid = [locX(1), locZ(1)];
    rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    % UM small anomalies
    rhopos(2).mid = rhoneg(1).mid + [-0.25 0.90] .* spacing;
    rhoneg(2).mid = rhoneg(1).mid + [+0.25 0.90] .* spacing;
    vspos(2).mid = vsneg(1).mid + [-0.25 0.90] .* spacing;
    vsneg(2).mid = vsneg(1).mid + [+0.25 0.90] .* spacing;
    vppos(2).mid = vpneg(1).mid + [-0.25 0.90] .* spacing;
    vpneg(2).mid = vpneg(1).mid + [+0.25 0.90] .* spacing;
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
            dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2) ]); 
            dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
            distUM_pos_rho = min([ norm([ii,jj] - rhopos(2).mid , 2) ]); 
            distUM_neg_rho = min([ norm([ii,jj] - rhoneg(2).mid , 2) ]); 
            distUM_pos_vs = min([ norm([ii,jj] - vspos(2).mid , 2) ]); 
            distUM_neg_vs = min([ norm([ii,jj] - vsneg(2).mid , 2) ]); 
            distUM_pos_vp = min([ norm([ii,jj] - vppos(2).mid , 2) ]); 
            distUM_neg_vp = min([ norm([ii,jj] - vpneg(2).mid , 2) ]); 
            
            % rho2 LM
            if dist_pos_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_rho < anom.radius;
                rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
            end
            % rho2 UM
            if distUM_pos_rho < anomUM.radius
                rho2(ii,jj) = rho2(ii,jj) * (1 + anomUM.strength);
            end
            if distUM_neg_rho < anomUM.radius
                rho2(ii,jj) = rho2(ii,jj) * (1 - anomUM.strength);
            end
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vs UM
            if distUM_pos_vs < anomUM.radius
                vs(ii,jj) = vs(ii,jj) * (1 + anomUM.strength);
            end
            if distUM_neg_vs < anomUM.radius
                vs(ii,jj) = vs(ii,jj) * (1 - anomUM.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
            % vp UM
            if distUM_pos_vp < anomUM.radius
                vp(ii,jj) = vp(ii,jj) * (1 + anomUM.strength);
            end
            if distUM_neg_vp < anomUM.radius
                vp(ii,jj) = vp(ii,jj) * (1 - anomUM.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2); 
    
    
elseif (model_type==86) % CIRCLES and small UM circles - ONLY VS AND VP:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
     
    % load prem
    [rho2, vs, vp] = load_PREM();
    
    % anomaly grid
    spacing = [nx / (3+2), nz / (2+1)]; % could change this to min(xspacing, zspacing) to have a regular grid
    %- location shift
    shift = [0.5 -0.2] .* spacing; % [to the right, up];
    %- divide into 3+1 horizontally
    locX = spacing(1) * [1:3] + shift(1);
    %- divide into 2+1 vertically
    locZ = spacing(2) * [1:2] + shift(2);
    [nx, nz];

    
    % anomaly properties
    anom.strength = 1.0 * 1/100; % 1 procent
    anom.radius = 0.30 * min(nx/(3+1), nz/(2+1));
%     anom.radius = 0.25 * min( locX(2)-locX(1), locZ(2)-locZ(1) ); % (2)-(1) so that I can still pad the outside of the domain
    anomUM.strength = anom.strength;
    anomUM.radius = 0.5 * anom.radius;
    

    % add the anomalies
    % --> you can probably do this without loops using any()
    %- rho at columns 1 and 4
%     rhopos(1).mid = [locX(1), locZ(1)];
%     rhoneg(1).mid = [locX(1), locZ(2)];
    %- vs at columns 2 & 5
    vspos(1).mid = [locX(2), locZ(1)];
    vsneg(1).mid = [locX(2), locZ(2)];
    %- vp at columns 3 & 6
    vppos(1).mid = [locX(3), locZ(1)];
    vpneg(1).mid = [locX(3), locZ(2)];
    
    % UM small anomalies
%     rhopos(2).mid = rhoneg(1).mid + [-0.25 0.90] .* spacing;
%     rhoneg(2).mid = rhoneg(1).mid + [+0.25 0.90] .* spacing;
    vspos(2).mid = vsneg(1).mid + [-0.25 0.90] .* spacing;
    vsneg(2).mid = vsneg(1).mid + [+0.25 0.90] .* spacing;
    vppos(2).mid = vpneg(1).mid + [-0.25 0.90] .* spacing;
    vpneg(2).mid = vpneg(1).mid + [+0.25 0.90] .* spacing;
    
    for ii = 1:size(rho2,1)
        for jj = 1:size(rho2,2)
%             dist_pos_rho = min([ norm([ii,jj] - rhopos(1).mid , 2) ]); 
%             dist_neg_rho = min([ norm([ii,jj] - rhoneg(1).mid , 2)]); 
            dist_pos_vs = min([ norm([ii,jj] - vspos(1).mid , 2)]); 
            dist_neg_vs = min([ norm([ii,jj] - vsneg(1).mid , 2)]); 
            dist_pos_vp = min([ norm([ii,jj] - vppos(1).mid , 2) ]); 
            dist_neg_vp = min([ norm([ii,jj] - vpneg(1).mid , 2)]); 
            
%             distUM_pos_rho = min([ norm([ii,jj] - rhopos(2).mid , 2) ]); 
%             distUM_neg_rho = min([ norm([ii,jj] - rhoneg(2).mid , 2) ]); 
            distUM_pos_vs = min([ norm([ii,jj] - vspos(2).mid , 2) ]); 
            distUM_neg_vs = min([ norm([ii,jj] - vsneg(2).mid , 2) ]); 
            distUM_pos_vp = min([ norm([ii,jj] - vppos(2).mid , 2) ]); 
            distUM_neg_vp = min([ norm([ii,jj] - vpneg(2).mid , 2) ]); 
            
%             % rho2 LM
%             if dist_pos_rho < anom.radius;
%                 rho2(ii,jj) = rho2(ii,jj) * (1 + anom.strength);
%             end
%             if dist_neg_rho < anom.radius;
%                 rho2(ii,jj) = rho2(ii,jj) * (1 - anom.strength);
%             end
%             % rho2 UM
%             if distUM_pos_rho < anomUM.radius
%                 rho2(ii,jj) = rho2(ii,jj) * (1 + anomUM.strength);
%             end
%             if distUM_neg_rho < anomUM.radius
%                 rho2(ii,jj) = rho2(ii,jj) * (1 - anomUM.strength);
%             end
            % vs
            if dist_pos_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vs < anom.radius;
                vs(ii,jj) = vs(ii,jj) * (1 - anom.strength);
            end
            % vs UM
            if distUM_pos_vs < anomUM.radius
                vs(ii,jj) = vs(ii,jj) * (1 + anomUM.strength);
            end
            if distUM_neg_vs < anomUM.radius
                vs(ii,jj) = vs(ii,jj) * (1 - anomUM.strength);
            end
            % vp
            if dist_pos_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 + anom.strength);
            end
            if dist_neg_vp < anom.radius;
                vp(ii,jj) = vp(ii,jj) * (1 - anom.strength);
            end
            % vp UM
            if distUM_pos_vp < anomUM.radius
                vp(ii,jj) = vp(ii,jj) * (1 + anomUM.strength);
            end
            if distUM_neg_vp < anomUM.radius
                vp(ii,jj) = vp(ii,jj) * (1 - anomUM.strength);
            end
        end
    end
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2); 
    
    
    
%% misc models    
elseif (model_type==100) % layered: left = high density, right = low density.
    
    rho=3000.0*ones(nx,nz);
    mu=4.8e10*ones(nx,nz);
    lambda=mu;
    
    rho(1:floor(nx/2),:) = 3500;
    
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
    
    
elseif (model_type==103) % PREM background model + raised 670 (by 30 km)
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!

    [rho2, vs, vp] = load_PREM();
    
    % how big is the deflection?
    input_parameters;
    [X,Z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
    width_defl = nx / 4;
    deflect670 = 30e3;
    low = min(deflect670, 0);
    high = max(deflect670,0);
        
    % do something to find the indices corresponding to the extent of the
    % 660 deflection
    left = floor(nx/2 - 0.5 * width_defl);
    right = ceil(nx/2 + 0.5 * width_defl);
    bips = Z(:,1);
    lo = find(bips >= (Lz - 670e3 + low));                 
    hi = find(bips <= (Lz - 670e3 + high));    
    bot = min(lo(1), hi(end));
    top = max(lo(1), hi(end));
    
    % do something to find the sub-660 values for rho2, vs, vp
    if deflect670 < 0
        s660rho = 3.992e3; s660vs = 5.570e3; s660vp = 10.266e3;
    else
        s660rho = 4.38074e3; s660vs = 5.94513e3; s660vp = 10.75132e3;
    end
    
    % adapt parameters within these boundaries to sub-660 values
    rho2(left:right, bot:top) = s660rho;
    vs(left:right, bot:top) = s660vs;
    vp(left:right, bot:top) = s660vp;
    
    % recalculating to rho-mu-lambda
    rho     = rho2;
    mu      = vs .^ 2 .* rho2;
    lambda  = rho2 .* ( vp.^2 - 2* vs.^2);
    
%% no model has been given    
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
    PREM = dlmread('./models/PREM-reference-model/PREM_1s_nowater.csv', '\t');
    
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

function mod_out = add_10randanoms(mod, anommax, whichrandmod)

nx = size(mod,1);
nz = size(mod,2);

% random locations of anomalies
    willekeurig(1).x = [0.906812, 0.493549, 0.499576, 0.175434, 0.441924, ...
                     0.668359, 0.541729, 0.444079, 0.348643, 0.927618];
    willekeurig(1).z = [0.261980, 0.930472, 0.645794, 0.534028, 0.907064, ...
                     0.631545, 0.039304, 0.194939, 0.657794, 0.749388];
    willekeurig(2).x = [0.814723686393179,0.905791937075619,...
                       0.126986816293506,0.913375856139019,...
                       0.632359246225410,0.0975404049994095,...
                       0.278498218867048,0.546881519204984,...
                       0.957506835434298,0.964888535199277];
    willekeurig(2).z = [0.157613081677548,0.970592781760616,...
                        0.957166948242946,0.485375648722841,...
                        0.800280468888800,0.141886338627215,...
                        0.421761282626275,0.915735525189067,...
                        0.792207329559554,0.959492426392903];
    willekeurig(3).x = [0.655740699156587,0.0357116785741896,...
                        0.849129305868777,0.933993247757551,...
                        0.678735154857774,0.757740130578333,...
                        0.743132468124916,0.392227019534168,...
                        0.655477890177557,0.171186687811562];
    willekeurig(3).z = [0.706046088019609,0.0318328463774207,...
                        0.276922984960890,0.0461713906311539,...
                        0.0971317812358475,0.823457828327293,...
                        0.694828622975817,0.317099480060861,...
                        0.950222048838355,0.0344460805029088];
                

wrm = whichrandmod;                 
    % add random anomalies to model
    for ii = 1:size(willekeurig(wrm).x,2)
        % anomaly strength
%         anom{i}.strength = 0.10 * 2600;
        anom{ii}.strength = anommax;
        
        anom{ii}.dxperc = willekeurig(wrm).x(ii);
        anom{ii}.dzperc = willekeurig(wrm).z(ii);
       
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