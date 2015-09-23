function Model_perc = make_model_percentage_of_mod86(percentageUM, percentageLM)
   
    
   % CIRCLES and small UM circles - ONLY VS AND VP:
                        % PREM background model plus regular grid of 'hard'
                        % edged circles - non-overlapping rho - vp - vs
                        % ONLY 6 CIRCLES
                        % IMPORTANT: 
                        % Model values will be sampled at height above CMB!
                        % so don't make the model higher than 2891 km!!
    
    % input
    input_parameters;
    
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
    
    % change anomaly strength using percentage
    anom.strength = percentageLM/100. * anom.strength;
    anomUM.strength = percentageUM/100. * anomUM.strength;
    

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
    Model_perc.rho     = rho2;
    Model_perc.mu      = vs .^ 2 .* rho2;
    Model_perc.lambda  = rho2 .* ( vp.^2 - 2* vs.^2); 
    
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