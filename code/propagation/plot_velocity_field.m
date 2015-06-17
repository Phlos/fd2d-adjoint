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

%%% Not happy with code structure -- too many ifs, not clear. NB 11-2-2014.

%%% Updated the way things are plotted -- no more ugly ifs 
%%%   --- Nienke Blom 18-4-2014.
    
path(path,'./tools')
figure(fig_vel)
clf

% Check what is calculated and what needs to be plotted
% if(strcmp(plot_forward_frames,'X-Y-Z'))
    
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
    
%- plot X, Z and Y velocity fields in a loop -----------------------------

for i=1:nplots
    
    subplot(1,nplots,i)
    hold on
    
    %- actual plotting -
    pcolor(X,Z,plotting(i).v');
    axis image
    
    %- colour scale -------------------------------------------------------
    cmin = -1*plotting(i).scale;
    cmax = plotting(i).scale;
    
    caxis([cmin cmax]);
    % flipud: colormap reverse to blue=down, red=up - it's not tomography..
    colormap(flipud(cm));
    shading interp
    
    
    %- axis labels and titles ---------------------------------------------
    
    xlabel('x [m]');
    ylabel('z [m]');
    
    title(plotting(i).title);
    
    if (strcmp(simulation_mode,'correlation') && t(n)<0)
        title('acausal correlation field');
    elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
        title('causal correlation field');
    end
    
    
    
    %======================================================================
    %- plot timestamp and source and receiver positions -------------------
    
    timestamp=['t [s] = ',num2str(floor(n*dt),'%4.0f')];
    timestamp_text = text(0.05*Lx,0.92*Lz,timestamp) ;
    text(0.95*Lx,0.92*Lz,['max = \pm', num2str(plotting(i).scale,'%3.1e')], ...
        'HorizontalAlignment','right')
    
    
    for p=1:5
        
        for k=1:length(src_x)
            plot(src_x(k),src_z(k),'kx')
        end
        
        for k=1:length(rec_x)
            plot(rec_x(k),rec_z(k),'ko')
        end
        
    end
    
    
    
end

hold off
    
    
    
    
% elseif(strcmp(plot_forward_frames,'PSV-SH'))
%     if(strcmp(wave_propagation_type,'SH'))
%         nplots=1;
%     elseif(strcmp(wave_propagation_type,'PSV'))
%         nplots=1;
%     elseif(strcmp(wave_propagation_type,'both'))
%         nplots=2;
%     end
%     
%     %- plot P-SV and SH velocity fields in a loop -----------------------------
%     
%     % P-SV velocity is just magnitude not direction (i.e. all positive)
%     vPSV = sqrt(vx.^2 + vz.^2);
%     
%     for i=1:nplots
%         subplot(1,nplots,i)
%         hold on
%         if(i==1)
%             if(strcmp(wave_propagation_type,'SH'))
%                 pcolor(X,Z,vy');
%             elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                 pcolor(X,Z,vPSV');
%             end
%         elseif(i==2)
%             pcolor(X,Z,vy');
%         end
%         axis image
%         
%         %- colour scale -------------------------------------------------------
%         
%         if(i==1)
%             if(strcmp(wave_propagation_type,'SH'))
%                 scale=prctile(abs(vy(:)),99.97);
%             elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                 scale=prctile(abs(vPSV(:)),99.97);
%             end
%         elseif(i==2)
%             scale=prctile(abs(vy(:)),99.97);
%         end
% 
%         
%         caxis([-scale scale]);
%         % flipud: colormap reverse to blue=down, red=up - it's not tomography..
%         colormap(flipud(cm));
%         shading interp
%         
%         
%         %- axis labels and titles ---------------------------------------------
%         
%         xlabel('x [m]');
%         ylabel('z [m]');
%         
%         if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))
%             if(i==1)
%                 if(strcmp(wave_propagation_type,'SH'))
%                     title('SH velocity field [m/s]');
%                 elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                     title('P-SV velocity field [m/s]');
%                 end
%             elseif(i==2)
%                 title('Y (SH) velocity field [m/s]');
%             end
%         elseif (strcmp(simulation_mode,'correlation') && t(n)<0)
%             title('acausal correlation field');
%         elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
%             title('causal correlation field');
%         end
%         
%         
%         
%         %======================================================================
%         %- plot timestamp and source and receiver positions -------------------
%         
%         timestamp=['t [s] = ',num2str(n*dt,'%4.0f')];
%         timestamp_text = text(0.05*Lx,0.92*Lz,timestamp) ;
%         text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ...
%             'HorizontalAlignment','right')
%         
%         
%         for p=1:5
%             
%             for k=1:length(src_x)
%                 plot(src_x(k),src_z(k),'kx')
%             end
%             
%             for k=1:length(rec_x)
%                 plot(rec_x(k),rec_z(k),'ko')
%             end
%             
%         end
%         
%         
%         
%     end
%     
%     hold off
%     
% elseif(strcmp(plot_forward_frames,'X-Y'))
%     
%     if(strcmp(wave_propagation_type,'SH'))
%         nplots=1;
%     elseif(strcmp(wave_propagation_type,'PSV'))
%         nplots=1;
%         vPSV=[vx vz];
%     elseif(strcmp(wave_propagation_type,'both'))
%         nplots=2;
%         vPSV=[vx vz];
%     end
%     
%     %- plot X, Z and Y velocity fields in a loop -----------------------------
%     
%     for i=1:nplots
%         subplot(1,nplots,i)
%         hold on
%         if(i==1)
%             if(strcmp(wave_propagation_type,'SH'))
%                 pcolor(X,Z,vy');
%             elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                 pcolor(X,Z,vx');
%             end
%         elseif(i==2)
%             pcolor(X,Z,vy');
%         end
%         axis image
%         
%         %- colour scale -------------------------------------------------------
%         
%         if(i==1)
%             if(strcmp(wave_propagation_type,'SH'))
%                 scale=prctile(abs(vy(:)),99.97);
%             elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                 scale=prctile(abs(vPSV(:)),99.97);
%             end
%         elseif(i==2)
%             scale=prctile(abs(vPSV(:)),99.97);
%         elseif(i==3)
%             scale=prctile(abs(vy(:)),99.97);
%         end
%         caxis([-scale scale]);
%         
%         % flipud: colormap reverse to blue=down, red=up - it's not tomography..
%         colormap(flipud(cm));
%         shading interp
%         
%         
%         %- axis labels and titles ---------------------------------------------
%         
%         xlabel('x [m]');
%         ylabel('z [m]');
%         
%         if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))
%             if(i==1)
%                 if(strcmp(wave_propagation_type,'SH'))
%                     title('Y (SH) velocity field [m/s]');
%                 elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%                     title('X velocity field [m/s]');
%                 end
%             elseif(i==2)
%                 title('Y (SH) velocity field [m/s]');
%             end
%         elseif (strcmp(simulation_mode,'correlation') && t(n)<0)
%             title('acausal correlation field');
%         elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
%             title('causal correlation field');
%         end
%         
%         
%         
%         %======================================================================
%         %- plot timestamp and source and receiver positions -------------------
%         
%         timestamp=['t [s] = ',num2str(n*dt,'%4.0f')];
%         timestamp_text = text(0.05*Lx,0.92*Lz,timestamp) ;
%         text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ...
%             'HorizontalAlignment','right')
%         
%         
%         for p=1:5
%             
%             for k=1:length(src_x)
%                 plot(src_x(k),src_z(k),'kx')
%             end
%             
%             for k=1:length(rec_x)
%                 plot(rec_x(k),rec_z(k),'ko')
%             end
%             
%         end
%         
%         
%         
%     end
    
    hold off
    
% end


%- record movie -------------------------------------------------------

if strcmp(make_movie,'yes')
    
    if exist('movie_index','var')
        movie_index=movie_index+1;
    else
        movie_index=1;
    end
    
    M(movie_index)=getframe(gcf);
    
end


%- save snapshot of wave propagation ----------------------------------

temps = n*dt;

if any(savetimes == temps)
    filename = [snapshotfile,'_forward_t',num2str(temps),'.png'];
    disp(['saving at time ',num2str(temps),' with filename ',filename, '!'])
    print(gcf,'-dpng','-r1000',filename)
end

%- finish -------------------------------------------------------------

pause(0.01)

