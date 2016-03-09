
function [m] = map_parameters_to_m(Model, usr_par)
% MAP_PARAMETERS_TO_M function to map the physical parameters rho, lambda
% and mu to the model variable m.
%
% Input: 
% rho 
% lambda
% mu
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% m
%
% See also MAP_M_TO_PARAMETERS and MAP_GRADPARAMETERS_TO_GRADM.

% input
input_parameters;

parametrisation = usr_par.parametrisation;
Model_bg        = usr_par.Model_bg;

% fixing vs and vp: there are only 1*nx*nz free parameters: rho
if strcmp(fix_velocities,'yes')
    
%    if strcmp(parametrisation, 'rhomulambda')
        m = Model.rho(:) ./ Model_bg.rho(:) - 1;
%    else
%        error('parametrisation must be rhomulambda if fixing velocities');
%    end

% no fixing of parameters: there are 3*nx*nz free parameters
else

    switch parametrisation
        
        % relative rho-mu-lambda (relative to bg model Model_bg)
        case 'rhomulambda'
            m1 = Model.rho ./ Model_bg.rho -1;
            m2 = Model.mu ./ Model_bg.mu -1;
            m3 = Model.lambda ./ Model_bg.lambda -1;
            
            m = [m1(:); m2(:); m3(:)];
            
        % relative rho-vs-vp (relative to bg model Model_bg)
        case 'rhovsvp'
            Mod_rvv = change_parametrisation('rhomulambda',parametrisation, Model);
            Mod_bg_rvv = change_parametrisation('rhomulambda',parametrisation, Model_bg);
            
            m4 = Mod_rvv.rho ./ Mod_bg_rvv.rho -1;
            m5 = Mod_rvv.vs ./ Mod_bg_rvv.vs -1;
            m6 = Mod_rvv.vp ./ Mod_bg_rvv.vp -1;
            
            m = [m4(:); m5(:); m6(:)];
            
    end
end

end
