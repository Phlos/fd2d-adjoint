function [outrho, fig_rhoupdate, drho, propsN ] = ...
                               apply_hard_constraints(propsR, inrho, axrot)

% function that calculates the new density model based on model properties
% of the current model and the known properties of the real model.
%
% INPUT:
% - propsR:         struct containing the properties of the real model:
%                   mass .M and moment of inertia .I (and model-insensitive
%                   properties such as volume .V, r^2vol .W and r^4vol .Z
% - inrho:          the test density model
%
% OUTPUT:
% - outrho:         output density model after applying hard constraints

%== PREPARATION ===========================================================

input_parameters;
[X,Zet,~,~]=define_computational_domain(Lx,Lz,nx,nz);

%- calculate model total moment of inertia I
if (axrot == 'x')
    r = abs(Zet');
elseif(axrot == 'z')
    r = abs(X');
else
    error('the rotation axis was not recognised')
end

%- calculate properties of current model
propsT = calculate_model_properties(inrho,axrot);

%- generalise W,V and Z
W = propsR.W;
V = propsR.V;
Z = propsR.Z;
if((propsT.V ~= V) || (propsT.W ~= W) || (propsT.Z ~= Z))
    error('The fixed domain properties (e.g. volume) appear not to be fixed...')
end


%== CALCULATION ===========================================================
%- calculate difference in mass and MoI between real and test models
Mdif = propsR.M - propsT.M;
Idif = propsR.I - propsT.I;
Mdifrel = Mdif / propsR.M;
Idifrel = Idif / propsR.I;


%- calculate factors A and B (for overview)
facA = (Z*V - W^2)^-1;
facB = W^2 / Z;


%- calculate mass contribution
%  constant
drho.Mconst = (facA * Z) * Mdif;
%  r^2 dependency
drho.Mvar = -(facA * W * r.^2) * Mdif;


%- calculate MoI contribution
%  constant
drho.Iconst = (-facA * W) * Idif;
%  r^2 dependency
drho.Ivar = (facA*facB + Z^-1) * r.^2 * Idif;


%- calculate rho update and new model properties
drho.total = drho.Mconst*ones(size(inrho)) + drho.Mvar ...
           + drho.Iconst*ones(size(inrho)) + drho.Ivar;

outrho = inrho + drho.total;

propsN = calculate_model_properties(outrho,axrot);


           
%== OUTPUT ================================================================
fprintf('Real, test mass  : %12.5e, %12.5e\n', propsR.M, propsT.M );
fprintf('Real, test MoI   : %12.5e, %12.5e\n', propsR.I, propsT.I );
fprintf('Real, test ratio : %12.5f, %12.5f\n\n', propsR.Ifactor, propsT.Ifactor );

fprintf('M diff  : %6.1e,   relative: %7.2e\n', Mdif, Mdifrel );
fprintf('I diff  : %6.1e,   relative: %7.2e\n\n', Idif, Idifrel );

fprintf('M const       : %10.1e\n', drho.Mconst );
fprintf('M var max, min: %10.1e %10.1e \n', ...
                       max(drho.Mvar(:)), min(drho.Mvar(:)) );
fprintf('I const       : %10.1e\n',drho.Iconst);
fprintf('I var max, min: %10.1e %10.1e \n\n', ...
                       max(drho.Ivar(:)), min(drho.Ivar(:)) );

fprintf('Max, min density update: %.2f, %.2f kg/m^3\n\n', ...
      max(max(outrho - inrho)), min(min(outrho - inrho)) );


fprintf('Real, new mass   : %12.5e, %12.5e\n', propsR.M, propsN.M );
fprintf('Real, new MoI    : %12.5e, %12.5e\n\n', propsR.I, propsN.I );

%- test whether the output density model gives us the real M and I
if((propsN.M - propsR.M > 3e-5) || (propsN.I - propsR.I > 3e-5))
    disp('ERRORR the updated mass and MoI are NOT the real mass and MoI!!');
    fprintf('Relative difference in mass, MoI: %10.2e, %10.2e \n\n', ...
              (propsR.M-propsN.M)/propsR.M, (propsR.I-propsN.I)/propsR.I );

end

fig_rhoupdate = figure;
load 'propagation/cm_model.mat';
colormap(cm_model);

subplot(2,2,1);
h = pcolor(X,Zet,(r.^2)');
shading interp
axis image;
set(h, 'EdgeColor', 'none');
caxis([-max(max(abs(r.^2))) max(max(abs(r.^2))) ]);
colorbar;
title('r^2');

subplot(2,2,2);
field = drho.Mvar + drho.Mconst*ones(size(inrho)) ;
i = pcolor(X,Zet,field');
shading interp
axis image;
set(i, 'EdgeColor', 'none');
cmax=abs(max(field(:)));
caxis([-cmax cmax]);
colorbar;
title('total M update');

subplot(2,2,3)
j = pcolor(X, Zet, (drho.Ivar + drho.Iconst*ones(size(inrho)) )' );
shading interp
axis image;
set(j, 'EdgeColor', 'none');
caxis([-max(max(abs(drho.Ivar))) max(max(abs(drho.Ivar))) ]);
colorbar;
title('total I update');

subplot(2,2,4);
k = pcolor(X, Zet, (drho.total)' );
shading interp
axis image;
set(k, 'EdgeColor', 'none');
caxis([-max(max(abs(drho.total))) max(max(abs(drho.total))) ]);
colorbar
title('total update')

end