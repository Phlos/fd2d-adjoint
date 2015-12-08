function properties = calc_mass_moi(rho, X, Z)
    % calculates the mass and inertia tensor for a given rho distribution
    
    %  0a) test if all fields are in the right shape
    % if (size(rho) ~= size(X)) || (size(rho) ~= size(Z))
    if ~( isequal(size(rho), size(X)) || isequal(size(rho), size(Z)) )
        X = X'; Z = Z';
    end
    if ~( isequal(size(rho), size(X)) || isequal(size(rho), size(Z)) )
        error('fields are not of the same size');
    end
    
    %  0b) determine dx, dz
    dx = X(2,1) - X(1,1);
    dz = Z(1,2) - Z(1,1);
    
    % R matrix
    R_11 = Z .^2;
    R_12 = -X .* Z;
    R_22 = X .^2;
    
    % calculate properties
    properties.mass = sum(dx*dz*rho(:));
    properties.I_11 = sum(R_11(:) .* rho(:) * dx * dz);
    properties.I_12 = sum(R_12(:) .* rho(:) * dx * dz);
    properties.I_22 = sum(R_22(:) .* rho(:) * dx * dz);
    
end