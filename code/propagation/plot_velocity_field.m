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

%%% Can do much better than this, all the stupid ifs. Instead of doing
%%% that, I could just make title objects and maybe also something for the
%%% nplots nonsense. Oh wellllll. not for now 11-2-2014.
    
    figure(fig_vel)
    clf
    
    if(strcmp(wave_propagation_type,'SH'))
        nplots=1;
    elseif(strcmp(wave_propagation_type,'PSV'))
        nplots=2;
    elseif(strcmp(wave_propagation_type,'both'))
        nplots=3;
    end
        
    %- plot X and Z velocity fields in a loop -----------------------------
    
    for i=1:nplots
    subplot(1,nplots,i)
    hold on
    if(i==1)
        if(strcmp(wave_propagation_type,'SH'))
            pcolor(X,Z,vy');
        elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
            pcolor(X,Z,vx');
        end
    elseif(i==2)
        pcolor(X,Z,vz');
    elseif(i==3)
        pcolor(X,Z,vy');
    end
    set(gca,'FontSize',20);
    axis image
    
    %- colour scale -------------------------------------------------------
    
%     if (n<0.8*length(t))
        if(i==1)
            if(strcmp(wave_propagation_type,'SH'))
                scale=max(max(abs(vy)));
            elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
                scale=max(max(abs(vx)));
            end
        elseif(i==2)
            scale=max(max(abs(vz)));
        elseif(i==3)
            scale=max(max(abs(vy)));
        end
%     end
   
    caxis([-scale scale]);
    % flipud: colormap reverse to blue=down, red=up - it's not tomography..
    colormap(flipud(cm));       
    shading interp

    
    %- axis labels and titles ---------------------------------------------
    
    xlabel('x [m]','FontSize',20);
    ylabel('z [m]','FontSize',20);
    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))
        if(i==1)
            if(strcmp(wave_propagation_type,'SH'))
                title('Y (SH) velocity field [m/s]','FontSize',20);
            elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
                title('X velocity field [m/s]','FontSize',20);
            end
        elseif(i==2)
            title('Z velocity field [m/s]','FontSize',20);
        elseif(i==3)
            title('Y (SH) velocity field [m/s]','FontSize',20);
        end
    elseif (strcmp(simulation_mode,'correlation') && t(n)<0)
        title('acausal correlation field','FontSize',20);
    elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
        title('causal correlation field','FontSize',20);
    end
    
    
    
    %======================================================================
    %- plot timestamp and source and receiver positions -------------------

        timestamp=['t [s] = ',num2str(n*dt)];
        timestamp_text = text(0.05*Lx,0.92*Lz,timestamp) ;
        text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ... 
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
    
%     %- save figure als jpeg -- takes just as much time as make_movie.
%     jpegname=['../output/jpeg/',num2str(n),'.jpg'];
%     print('-djpeg', jpegname);
    
    %- record movie -------------------------------------------------------
    
    if strcmp(make_movie,'yes')
    
        if exist('movie_index','var')
            movie_index=movie_index+1;
        else
            movie_index=1;
        end
        
        M(movie_index)=getframe(gcf);
        
    end
    
    %- finish -------------------------------------------------------------
    
    pause(0.01)
    
