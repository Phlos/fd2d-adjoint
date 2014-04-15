% calculate other kernels from the basis kernels in rho mu lambda

function [K, K_rel] = calculate_other_kernels(K)

path(path,'../input')
input_parameters;

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
[mu,rho,lambda]=define_material_parameters(nx,nz,11);

%% fill out the kernels which are not calculated but which one may want to plot

if (strcmp(wave_propagation_type,'SH'))
    K.mu.PSV = zeros(size(K.mu.SH));
    K.rho.PSV = zeros(size(K.rho.SH));
    K.lambda.PSV = zeros(size(K.rho.SH));
end

if (strcmp(wave_propagation_type,'PSV'))
    K.mu.SH = zeros(size(K.mu.PSV));
    K.rho.SH = zeros(size(K.rho.PSV));
end

%- total kernels, used for the inversion based on rho mu lambda
%  parametrisation.
K.rho.total = K.rho.PSV + K.rho.SH;
K.mu.total = K.mu.PSV + K.mu.SH;
K.lambda.total = K.lambda.PSV;

%% calculate kappa, vp, vs for the model
kappa = (lambda + 2/3 * mu);
vs = sqrt( mu ./ rho );
vp = sqrt( (lambda + 2*mu) ./ rho );

%% parametrisation rho, mu, kappa
K.rho1.PSV = K.rho.PSV;
K.mu1.PSV = K.mu.PSV - 2/3* K.lambda.PSV;
K.kappa1.PSV = K.lambda.PSV;

K.rho1.SH = K.rho.SH;
K.mu1.SH = K.mu.SH; % - 2/3* K.lambda.SH;
% K.kappa1.SH = K.lambda.SH;                  % should be zero anyway!!

%% parametrisation rho, vs. vp



% actual kernels
K.rho2.PSV = K.rho.PSV + (vp.^2 - 2.*vs.^2) .* K.lambda.PSV ...
                       + vs.^2 .* K.mu.PSV;
K.vs2.PSV = 2*rho .* vs .* K.mu.PSV - 4 * rho .* vs .* K.lambda.PSV;
K.vp2.PSV = 2*rho .* vp .* K.lambda.PSV;

K.rho2.SH = K.rho.SH ... % + (vp.^2 - 2.*vs.^2) .* K.lambda.SH ...
                       + vs.^2 .* K.mu.SH;
K.vs2.SH = 2*rho .* vs .* K.mu.SH; %- 4 * rho .* vs .* K.lambda.SH;
% K.vp2.SH = 2*rho .* vp .* K.lambda.SH;



%% relative kernels
K_rel.rho.PSV = K.rho.PSV.*rho;
K_rel.rho.SH = K.rho.SH.*rho;
K_rel.mu.PSV = K.mu.PSV.*mu;
K_rel.mu.SH = K.mu.SH.*mu;
K_rel.lambda.PSV = K.lambda.PSV.*lambda;

K_rel.rho1.PSV = K.rho1.PSV.*rho;
K_rel.rho1.SH = K.rho1.SH.*rho;
K_rel.mu1.PSV = K.mu1.PSV.*mu;
K_rel.mu1.SH = K.mu1.SH.*mu;
K_rel.kappa1.PSV = K.kappa1.PSV .* kappa;

K_rel.rho2.PSV = K.rho2.PSV.*rho;
K_rel.rho2.SH = K.rho2.SH.*rho;
K_rel.vs2.PSV = K.vs2.PSV.*vs;
K_rel.vs2.SH = K.vs2.SH.*vs;
K_rel.vp2.PSV = K.vp2.PSV.*vp;

% 
% no_SH = {'lambda', 'kappa1', 'vp2'};
% fn_params = fieldnames(K);
% fn_dirs = {'PSV'; 'SH'};
% for i = 1:length(fn_params)
% %     fn_params{i}
%     for j = 1:length(fn_dirs)
%         disp([fn_params{i},'.',fn_dirs{j}])
%         bops = strcmp(fn_params{i},no_SH); % if one of these is 1 then don't calculate SH if dir is SH
%         niks = false(size(bops));
%         if ( strcmp(fn_dirs{j},'SH' )  )
%         Krel.(fn_params{i}).(fn_dirs{j}) = ...
%                  K.(fn_params{i}).(fn_dirs{j}) .* eval(fn_params{i});
%         else
%             disp(['we may not calculate K.',fn_params{i},'.',fn_dirs{j}])
%         end
%     end;
%     
% end