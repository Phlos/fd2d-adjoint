%==========================================================================
% Script to plot the velocity field of P-SV and/or SH in 2D
%
%   "INPUT": (meaning: variables that need to be present and useable)
% ---------
% - n:                      the timestep we're at.
% - wave_propagation_type:  'SH', 'PSV', 'both'
% - X, Z (of the plot):      x and z coordinates over the whole plot
% - vy, vx, vz:               the velocity fields. (not all needed -- 
%                           depends on wave_propagation_type)
% - simulation_mode:        'forward', 'forward_green', 'correlation', ...
% - make_movie:             'yes' or 'no'. If plotting only one - say no.
%
%   OUTPUT:
%-----------
% a figure with the requested velocity fields.
%==========================================================================

%%% Updated the way things are plotted -- no more ugly ifs 
%%%   --- Nienke Blom 18-4-2014.
    
% path(path,'./tools')
input_parameters;
figure(fig_vel)
clf

%% Plotting prep - check what is calculated and what needs to be plotted

    
if(strcmp(wave_propagation_type,'SH'))
    nplots=1;
    plotting(1).v       = vy;
    plotting(1).scale   = prctile(abs(vy(:)),99.97);
    plotting(1).title   = 'Y velocity field [m/s]';
    
    if strcmp(plot_forward_frames,'X-Y-Z')
        warning('Plotting:ignoringPSV', ...
            'Not plotting the X and Z velocity fields as they''re not calculated')
        
    elseif strcmp(plot_forward_frames,'X-Y')
        warning('Plotting:ignoringPSV', ...
            'Not plotting the X velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV-SH')
        plotting(1).title = 'SH velocity field [m/s]';
        warning('Plotting:ignoringPSV', ...
            'Not plotting the P-SV velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV')
        error('You cannot calculate SH and plot P-SV wave propagation')
    
    elseif strcmp(plot_forward_frames,'SH')
        plotting(1).title   = 'SH velocity field [m/s]';
    
    end
    
elseif(strcmp(wave_propagation_type,'PSV'))
    
    if strcmp(plot_forward_frames,'X-Y-Z')
        nplots=2;
        vxz = [vx vz];
        plotting(1).v       = vx;
        plotting(1).scale   = prctile(abs(vxz(:)),99.97);
        plotting(1).title   = 'X velocity field [m/s]';
        
        plotting(2).v       = vz;
        plotting(2).scale   = prctile(abs(vxz(:)),99.97);
        plotting(2).title   = 'Z velocity field [m/s]';
        warning('Plotting:ignoringSH', ...
            'Not plotting the Y velocity field as it''s not calculated')
        
    elseif strcmp(plot_forward_frames,'X-Y')
        nplots=1;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).scale   = prctile(abs(vxz(:)),99.97);
        plotting(1).title   = 'X velocity field [m/s]';
        warning('Plotting:ignoringSH', ...
            'Not plotting the Y velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV-SH')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        
        plotting(1).v       = vPSV;
        plotting(1).scale   = prctile(abs(vPSV(:)),99.97);
        plotting(1).title   = 'P-SV velocity field [m/s]';
        warning('Plotting:ignoringY', ...
            'Not plotting the SH velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        
        plotting(1).v       = vPSV;
        plotting(1).scale   = prctile(abs(vPSV(:)),99.97);
        plotting(1).title   = 'P-SV velocity field [m/s]';
    
    elseif strcmp(plot_forward_frames,'SH')
        error('You cannot calculate PSV and plot SH wave propagation')
    end
    
elseif(strcmp(wave_propagation_type,'both'))
    if strcmp(plot_forward_frames,'X-Y-Z')
        nplots=3;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).scale   = prctile(abs(vxz(:)),99.97);
        plotting(1).title   = 'X velocity field [m/s]';
        
        plotting(2).v       = vz;
        plotting(2).scale   = prctile(abs(vxz(:)),99.97);
        plotting(2).title   = 'Z velocity field [m/s]';

        plotting(3).v       = vy;
        plotting(3).scale   = prctile(abs(vy(:)),99.97);
        plotting(3).title   = 'Y velocity field [m/s]';

    elseif strcmp(plot_forward_frames,'X-Y')
        nplots=2;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).scale   = prctile(abs(vxz(:)),99.97);
        plotting(1).title   = 'X velocity field [m/s]';
        
        plotting(2).v       = vy;
        plotting(2).scale   = prctile(abs(vy(:)),99.97);
        plotting(2).title   = 'Y velocity field [m/s]';
        
    elseif strcmp(plot_forward_frames,'PSV-SH')
        nplots=2;
        vPSV = sqrt(vx.^2 + vz.^2);
        
        plotting(1).v       = vPSV;
        plotting(1).scale   = prctile(abs(vPSV(:)),99.97);
        plotting(1).title   = 'P-SV velocity field [m/s]';
        
        plotting(2).v       = vy;
        plotting(2).scale   = prctile(abs(vy(:)),99.97);
        plotting(2).title   = 'SH velocity field [m/s]';
        
    elseif strcmp(plot_forward_frames,'PSV')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        
        plotting(1).v       = vPSV;
        plotting(1).scale   = prctile(abs(vPSV(:)),99.97);
        plotting(1).title   = 'P-SV velocity field [m/s]';
        
    elseif strcmp(plot_forward_frames,'SH')
        nplots = 1;
        
        plotting(1).v       = vy;
        plotting(1).scale   = prctile(abs(vy(:)),99.97);
        plotting(1).title   = 'SH velocity field [m/s]';
    end
    
end
 

% convert distances to km & set Z to depth below surface

if any( model_type == [50 : 1 : 89] )
    surface_level = 2890; % [km] only valid in PREM
else
    surface_level = Lz / 1000;
end
X_new = X ./ 1000;
Z_new = Z ./ 1000;
Z_new = surface_level - (Z_new);
for k=1:length(src_info)
    src_x_new(k) = src_info(k).loc_x / 1000;
    src_z_new(k) = src_info(k).loc_z / 1000;
    src_z_new(k) = surface_level - src_z_new(k);
end

for k=1:length(rec_x)
    rec_x_new(k) = rec_x(k) ./ 1000;
    rec_z_new(k) = rec_z(k) ./ 1000;
    rec_z_new(k) = surface_level - rec_z_new(k);
end

%% - plot X, Z and Y velocity fields in a loop -----------------------------

for i=1:nplots
    
    subplot(1,nplots,i)
    hold on
    
    % make sure that the Y axis points downwards
    set(gca, 'YDir', 'reverse');
    
    %- actual plotting -
    pcolor(X_new,Z_new,plotting(i).v');
    axis image
    
    %- colour scale -------------------------------------------------------
    cmin = -1*plotting(i).scale;
    cmax = plotting(i).scale;
    
    caxis([cmin cmax]);
    % flipud: colormap reverse to blue=down, red=up - it's not tomography..
    colormap(flipud(cm));
    shading interp
    
    
    %- axis labels and titles ---------------------------------------------
    
    xlabel('x [km]');
    ylabel('z [km]');
    
    title(plotting(i).title);
    
    if (strcmp(simulation_mode,'correlation') && t(n)<0)
        title('acausal correlation field');
    elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
        title('causal correlation field');
    end
    
    
    
    %======================================================================
    %- plot source and receiver positions ---------------------------------
    
        
    for k=1:length(src_x_new)
        plot(src_x_new(k),src_z_new(k),'kx')
    end
    
    for k=1:length(rec_x_new)
        plot(rec_x_new(k),rec_z_new(k),'ko')
    end
    
    %- plot contour of rho anomaly
    if exist('plot_contour_rho_anom', 'var') && plot_contour_rho_anom
        contour(X_new, Z_new, model.rho', 1, '--k');
    end
    
    %- plot timestamp and velocity max value
    timestamp_text=['t = ',num2str(floor(n*dt),'%4.0f'), ' s'];
    timestamp_plot = text(0.05,0.92, timestamp_text, ...
        'HorizontalAlignment','left', 'Units', 'normalized');
%     text(0.05*Lx,0.92*Lz,timestamp_text) ;
    max_text = ['max = \pm', num2str(plotting(i).scale,'%3.1e')];
    text(0.95,0.92, max_text, ...
        'HorizontalAlignment','right', 'Units', 'normalized');
    
    %- plot wavefield description
    if exist('movie_label', 'var')
    text(0.95,0.08, movie_label, ...
        'HorizontalAlignment','right', 'Units', 'normalized');
    end
    
end

hold off


%% - record movie -------------------------------------------------------

if strcmp(make_movie,'yes')
    

    
    if exist('movie_index','var')
        movie_index=movie_index+1;
    else
        movie_index=1;
    end
    
    M(movie_index)=getframe(gcf);
    
end


%% - save snapshot of wave propagation ----------------------------------

temps = n*dt;

if any(savetimes == temps)
    filename = [snapshotfile,'_forward_t',num2str(temps),'.png'];
    disp(['saving at time ',num2str(temps),' with filename ',filename, '!'])
    print(gcf,'-dpng','-r1000',filename)
end

%- finish -------------------------------------------------------------

pause(0.01)

