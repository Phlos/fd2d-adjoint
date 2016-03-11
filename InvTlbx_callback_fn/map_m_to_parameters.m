
function [Model] = map_m_to_parameters(m, usr_par)
% MAP_M_TO_PARAMETERS function to map the model variable m to physical
% parameters rho, lambda and mu.
%
% Input: 
% m
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% rho 
% lambda
% mu
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.

% input
input_parameters;

parametrisation = usr_par.parametrisation;
Model_bg        = usr_par.Model_bg;
Model_start     = usr_par.Model_start;


% fixing vs and vp: there are only 1*nx*nz free parameters: rho
if strcmp(fix_velocities,'yes')
%     vs_bg = sqrt(Model_bg.mu ./ Model_bg.rho);
%     vp_bg = sqrt( (Model_bg.lambda + 2*Model_bg.mu) ./ Model_bg.rho);
    vs_st = sqrt(Model_start.mu ./ Model_start.rho);
    vp_st = sqrt( (Model_start.lambda + 2*Model_start.mu) ./ Model_start.rho);
    
%    if strcmp(parametrisation, 'rhomulambda')
        rho = Model_bg.rho .* (1 + reshape(m, nx, nz));
        mu  = vs_st .^2 .* rho;
        lambda = (vp_st.^2 - 2* vs_st.^2) .* rho;
            
%    else
%        error('parametrisation must be rhomulambda if fixing velocities');
%    end

% fixing density: there are 2*nx*nz free parameters: mu+lambda / vs+vp
elseif strcmp(fix_density, 'yes');

	switch parametrisation
		case 'rhomulambda';

			m2 = m(        1 :   nx*nz);
			m3 = m(nx*nz + 1 : 2*nx*nz);

			rho    = Model_start.rho;
			mu     = Model_bg.mu     .* (1+ reshape(m2, nx,nz));
			lambda = Model_bg.lambda .* (1+ reshape(m3, nx,nz));
		
		case 'rhovsvp';

			Model_bg = change_parametrisation('rhomulambda', 'rhovsvp', Model_bg);

			m5 = m(        1 :   nx*nz);
			m6 = m(nx*nz + 1 : 2*nx*nz);

			rho = Model_start.rho;
			vs  = Model_bg.vs .* (1+ reshape(m5, nx,nz));
			vp  = Model_bg.vp .* (1+ reshape(m6, nx,nz));

			mu = rho .* vs.^2;
			lambda = rho .* (vp .^2 - 2.*vs .^2);
	end


% no fixing of parameters: there are 3*nx*nz free parameters
else
    
    
    % determine which part of the m vector is which parameters
    switch parametrisation
        case 'rhomulambda'
            m1 = m(          1 :   nx*nz);
            m2 = m(  nx*nz + 1 : 2*nx*nz);
            m3 = m(2*nx*nz + 1 : 3*nx*nz);
            
            % reparametrise
            rho = Model_bg.rho .* (1 + reshape(m1, nx, nz));
            mu  = Model_bg.mu  .* (1 + reshape(m2, nx, nz));
            lambda = Model_bg.lambda .* (1 + reshape(m3, nx, nz));
            
        case 'rhovsvp'
            
            Model_bg = change_parametrisation('rhomulambda', 'rhovsvp', Model_bg);
            
            m4 = m(          1 :   nx*nz);
            m5 = m(  nx*nz + 1 : 2*nx*nz);
            m6 = m(2*nx*nz + 1 : 3*nx*nz);
            
            % reparametrise
            rho = Model_bg.rho .* (1 + reshape(m4, nx, nz));
            vs = Model_bg.vs  .* (1 + reshape(m5, nx, nz));
            vp = Model_bg.vp .* (1 + reshape(m6, nx, nz));
            
            mu = rho .* vs.^2;
            lambda = rho .* (vp .^2 - 2.*vs .^2);
            
    end
end

% reshape
Model.rho = reshape(rho,nx,nz);
Model.mu = reshape(mu,nx,nz);
Model.lambda = reshape(lambda,nx,nz);


end
