
function [gradm] = map_gradparameters_to_gradm(K_total, usr_par)
% MAP_GRADPARAMETERS_TO_GRADM function to transform the partial derivatives
% w.r.t. the physical parameters rho, lambda to the partial derivatives
% w.r.t. to the model variable m. This requires to apply the chain rule
% which is implicitly given by the function MAP_M_TO_PARAMETERS.
%
% Input: 
% gradrho 
% gradlambda
% gradmu
%
% Output:
% gradm
%
% See also MAP_M_TO_PARAMETERS.

%% input
% input
input_parameters;
Model = usr_par.Mod_current;
Model_bg = usr_par.Model_bg;

%% preparation of kernel
% reparametrising Kernel to inversion parametrisation
K_reparam = change_parametrisation_kernels('rhomulambda', parametrisation, K_total, Model);

% calculate relative kernels
K_rel = calculate_relative_kernels(K_reparam, Model_bg);

% set bottom five rows to zero if requested (5 rows = width FD stencil
% ( this is to avoid irritating focusing at the bottom of the domain due to
%   numerical instabilities of the FD scheme at boundary conditions )
n0rows = 1+5;
if strcmp(zero_bottom_rows, 'yeszerobottom')
    % set bottom rows to zero
    params = fieldnames(K_rel);
    for ii = 1:numel(params)
        comps = fieldnames(K_rel.(params{ii}));
        for jj = 1:numel(comps)
            K_rel.(params{ii}).(comps{jj})(:,1:n0rows) = 0;
        end
    end
end

% filter kernels
if strcmp(smoothing, 'nosmooth')
    disp('WARNING! Kernels are not being filtered!!');
else
    K_rel = filter_kernels(K_rel, parametrisation, smoothgwid);
end

%% actual mapping into inversion parametrisation
% fixing vs and vp: there are only 1*nx*nz free parameters: rho
if strcmp(fix_velocities,'yes')
    
    if strcmp(parametrisation, 'rhomulambda')
        gradm = [K_rel.rho.total(:) ...
            + K_rel.mu.total(:) ...
            + K_rel.lambda.total(:) ];
    else
	gradm = K_rel.rho2.total(:);
%        error('parametrisation must be rhomulambda if fixing velocities');
    end
    
% no fixing of parameters: there are 3*nx*nz free parameters
else

    switch parametrisation
        case 'rhomulambda'
            gradm = [K_rel.rho.total(:); ...
                K_rel.mu.total(:); ...
                K_rel.lambda.total(:)];
        case 'rhovsvp'
%            if strcmp(fix_velocities,'yes')
%                error('WARNING! no fix vp vs defined here!')
%            end
            gradm = [K_rel.rho2.total(:); ...
                K_rel.vs2.total(:); ...
                K_rel.vp2.total(:)];
    end


end

end
