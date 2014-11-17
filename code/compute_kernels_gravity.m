function [Kg, fig_Kg] = compute_kernels_gravity(g_src, rec_grav, varargin)

% calculation of gravity kernels based on the gravity source (obtained in
% e.g. calculate_gravity_misfit) and on domain properties.

%- PREPARATION

%- varargin
plotornot = checkargs(varargin(:)); % plots figures by default

%- prepare necessary information
path(path,'../input');
path(path,'../tools');
path(path,'../code');
path(path,'../code/propagation');
set_figure_properties_doffer;

% input
input_parameters;
[X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);


% gravitational constant
% 6.67384 × 10-11 m3 kg-1 s-2
G = 6.67384e-11;

nrec = size(rec_grav.x,2);


%- COMPUTATION


if(strcmp(plotornot,'yes'))
    fig_grav_src = figure;
    recs = [1:nrec];
    plot(recs, g_src.x, recs, g_src.z);
    %     close(gcf)
end

%- actual gravity kernel calculation

for i = 1:nrec % loop over recorders
    
    % calculate distance vector r{i}.x, r{i}.z
    r{i}.x = rec_grav.x(i) - X';
    r{i}.z = rec_grav.z(i) - Z';
    
    % calculate length or the vectors r for each point
    r{i}.length = sqrt(r{i}.x.^2 + r{i}.z.^2);
    
    % gravity kernel per component per receiver
    Kg_rec{i}.x = -G * g_src.x(i) * r{i}.x ./ r{i}.length .^ 3;
    Kg_rec{i}.z = -G * g_src.z(i) * r{i}.z ./ r{i}.length .^ 3;
    
    % gravity kernel per receiver
    Kg_rec{i}.total = Kg_rec{i}.x + Kg_rec{i}.z;
    
end



%- OUTPUT

% initialisation
Kg = zeros(size(Kg_rec{1}.total));

% initialise plot
if(strcmp(plotornot,'yes'))
load '../code/propagation/cm_velocity.mat';
    totalkernel = figure;
    set(totalkernel,'OuterPosition',pos_gravknl_buildup)
    % receiverkernels = figure;
end

for i = 1:nrec
    % kernel calculation
    Kg = Kg + Kg_rec{i}.total;
        
    if(strcmp(plotornot,'yes'))
        
        % plot total kernel buildup
        % top left
        figure(totalkernel)
        subplot(2,2,1)
        Kg_sm = filter_kernels(Kg, 15);
        pcolor(X,Z, Kg_sm')
        colormap(cm);
        shading interp
        axis image
        titel = ['total kernel'];
        title(titel);
        colorbar
        scale = prctile(abs(Kg_sm(:)),98);
        caxis([-scale scale]);
        % top right
        subplot(2,2,2)
%         Kg_rec_sm.total = filter_kernels(Kg_rec{i}.total, 15);
        pcolor(X,Z, Kg_rec{i}.total')
        shading interp
        axis image
        titel = ['total kernel (receiver ', num2str(i),')'];
        title(titel);
        colorbar
        scale = prctile(abs(Kg_rec{i}.total(:)),98);
        caxis([-scale scale]);
        % bottom left
        subplot(2,2,3)
%         Kg_rec_sm.x = filter_kernels(Kg_rec{i}.x, 15);
        pcolor(X,Z, Kg_rec{i}.x')
        shading interp
        axis image
        titel = ['x kernel (receiver ', num2str(i),')'];
        title(titel);
        colorbar
%         scale = prctile(abs(Kg_rec_sm.x(:)),98);
        caxis([-scale scale]);
        % bottom right
        subplot(2,2,4)
%         Kg_rec_sm.z = filter_kernels(Kg_rec{i}.z, 15);
        pcolor(X,Z, Kg_rec{i}.z')
        shading interp
        axis image
        titel = ['z kernel (receiver ', num2str(i),')'];
        title(titel);
        colorbar
%         scale = prctile(abs(Kg_rec_sm.z(:)),98);
        caxis([-scale scale]);
    
%     % plot of sub kernels and total kernel
%     figure(receiverkernels)
%     subplot(4, ceil(nrec/4),i)
%     pcolor(X,Z, Kg_rec{i}.total')
%     shading interp
%     axis image
%     colorbar
%     titel = ['total kernel (receiver ', num2str(i),')'];
%     title(titel);
    
%     pause(1);

    end
    
end

if(strcmp(plotornot,'yes'))
    close(fig_grav_src, totalkernel);
else
        Kg_sm = filter_kernels(Kg, 15); % if it hasn't been smoothed already, smooth now for output figure
end

%- output figure
fig_Kg = figure;
pcolor(X,Z, Kg_sm')
colormap(cm);
shading interp
axis image
titel = ['Total kernel'];
title(titel);
colorbar
scale = prctile(abs(Kg_sm(:)),98);
caxis([-scale scale]);
      

end

function plotornot = checkargs(arg)

% size(arg)
narg = size(arg,1);

if narg == 1
    plotornot = arg{1};
elseif narg == 0
    plotornot = 'yes';
else
    error('we don''t know whether you want to plot the g kernel or not')
end

end