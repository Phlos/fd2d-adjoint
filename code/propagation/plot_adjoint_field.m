% plot the adjoint fields
%==========================================================================

figure(fig_adjoint);
ncols=4;
% set color map
colormap(flipud(cm));

temps = nt*dt - (n-plot_every)*dt;
timestamp=['t [s] = ',num2str(temps,'%4.0f')];
% disp(timestamp);


%% Check what is calculated and what needs to be plotted
% if(strcmp(plot_forward_frames,'X-Y-Z'))
    
if(strcmp(wave_propagation_type,'SH'))
    nplots=1;
    plotting(1).v       = vy;
    plotting(1).v_fw    = vy_fw;
    plotting(1).interact= interaction.rho.y;
    plotting(1).kernel  = K.rho.total;
    plotting(1).dir     = 'SH';
    plotting(1).knlname = '\rho';
    
    if strcmp(plot_forward_frames,'X-Y-Z')
        warning('Plotting:ignoringPSV', ...
            'Not plotting the X and Z velocity fields as they''re not calculated')
        
    elseif strcmp(plot_forward_frames,'X-Y')
        warning('Plotting:ignoringPSV', ...
            'Not plotting the X velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV-SH')
%         plotting(1).title = 'SH velocity field [m/s]';
        warning('Plotting:ignoringPSV', ...
            'Not plotting the P-SV velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV')
        error('You cannot calculate SH and plot P-SV wave propagation')
    
    elseif strcmp(plot_forward_frames,'SH')
%         plotting(1).title   = 'SH velocity field [m/s]';
    
    end
    
elseif(strcmp(wave_propagation_type,'PSV'))
    
    if strcmp(plot_forward_frames,'X-Y-Z')
        nplots=2;
        vxz = [vx vz];
        plotting(1).v       = vx;
        plotting(1).v_fw    = vx_fw;
        plotting(1).interact= interaction.rho.x;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'X';
        plotting(1).knlname = '\rho';
        
        plotting(2).v       = vz;
        plotting(2).v_fw    = vz_fw;
        plotting(2).interact= interaction.rho.z;
        plotting(2).kernel  = K.mu.total;
        plotting(2).dir     = 'Z';
        plotting(2).knlname = '\mu';
        warning('Plotting:ignoringSH', ...
            'Not plotting the Y velocity field as it''s not calculated')
        
    elseif strcmp(plot_forward_frames,'X-Y')
        nplots=1;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).v_fw    = vx_fw;
        plotting(1).interact= interaction.rho.x;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'X';
        plotting(1).knlname = '\rho';
        warning('Plotting:ignoringSH', ...
            'Not plotting the Y velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV-SH')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        vPSV_fw = sqrt(vx_fw.^2 + vz_fw.^2);
        interaction_PSV = sqrt(interaction.rho.x .^2 + interaction.rho.z .^2);
        
        plotting(1).v       = vPSV;
        plotting(1).v_fw    = vPSV_fw;
        plotting(1).interact= interaction_PSV;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'P-SV';
        plotting(1).knlname = '\rho';
        warning('Plotting:ignoringY', ...
            'Not plotting the SH velocity field as it''s not calculated')
    
    elseif strcmp(plot_forward_frames,'PSV')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        vPSV_fw = sqrt(vx_fw.^2 + vz_fw.^2);
        interaction_PSV = sqrt(interaction.rho.x .^2 + interaction.rho.z .^2);
        
        plotting(1).v       = vPSV;
        plotting(1).v_fw    = vPSV_fw;
        plotting(1).interact= interaction_PSV;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'P-SV';
        plotting(1).knlname = '\rho';
    
    elseif strcmp(plot_forward_frames,'SH')
        error('You cannot calculate PSV and plot SH wave propagation')
    end
    
elseif(strcmp(wave_propagation_type,'both'))
    if strcmp(plot_forward_frames,'X-Y-Z')
        nplots=3;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).v_fw    = vx_fw;
        plotting(1).interact= interaction.rho.x;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'X';
        plotting(1).knlname = '\rho';
        
        plotting(2).v       = vz;
        plotting(2).v_fw    = vz_fw;
        plotting(2).interact= interaction.rho.z;
        plotting(2).kernel  = K.mu.total;
        plotting(2).dir     = 'Z';
        plotting(2).knlname = '\mu';
        plotting(3).v       = vy;
        plotting(3).v_fw    = vy_fw;
        plotting(3).interact= interaction.rho.y;
        plotting(3).kernel  = K.lambda.total;
        plotting(3).dir     = 'Y';
        plotting(3).knlname = '\lambda';

    elseif strcmp(plot_forward_frames,'X-Y')
        nplots=2;
        vxz = [vx vz];
        
        plotting(1).v       = vx;
        plotting(1).v_fw    = vx_fw;
        plotting(1).interact= interaction.rho.x;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'X';
        plotting(1).knlname = '\rho';
        
        plotting(2).v       = vy;
        plotting(2).v_fw    = vy_fw;
        plotting(2).interact= interaction.rho.y;
        plotting(2).kernel  = K.mu.total;
        plotting(2).dir     = 'Y';
        plotting(2).knlname = '\mu';
        
    elseif strcmp(plot_forward_frames,'PSV-SH')
        nplots=2;
        vPSV = sqrt(vx.^2 + vz.^2);
        vPSV_fw = sqrt(vx_fw.^2 + vz_fw.^2);
        interaction_PSV = sqrt(interaction.rho.x .^2 + interaction.rho.z .^2);
        
        plotting(1).v       = vPSV;
        plotting(1).v_fw    = vPSV_fw;
        plotting(1).interact= interaction_PSV;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'P-SV';
        plotting(1).knlname = '\rho';
        
        plotting(2).v       = vy;
        plotting(2).v_fw    = vy_fw;
        plotting(2).interact= interaction.rho.y;
        plotting(2).kernel  = K.mu.total;
        plotting(2).dir     = 'SH';
        plotting(2).knlname = '\mu';
        
    elseif strcmp(plot_forward_frames,'PSV')
        nplots=1;
        vPSV = sqrt(vx.^2 + vz.^2);
        vPSV_fw = sqrt(vx_fw.^2 + vz_fw.^2);
        interaction_PSV = sqrt(interaction.rho.x .^2 + interaction.rho.z .^2);
        
        plotting(1).v       = vPSV;
        plotting(1).v_fw    = vPSV_fw;
        plotting(1).interact= interaction_PSV;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'P-SV';
        plotting(1).knlname = '\rho';
        
    elseif strcmp(plot_forward_frames,'SH')
        nplots = 1;
        
        plotting(1).v       = vy;
        plotting(1).v_fw    = vy_fw;
        plotting(1).interact= interaction.rho.y;
        plotting(1).kernel  = K.rho.total;
        plotting(1).dir     = 'SH';
        plotting(1).knlname = '\rho';
    end
    
end


%% every row is a 'direction': x, y or z, P-SV or SH...
for i=1:nplots
%     if(row==1)
%         if(strcmp(wave_propagation_type,'SH'))
%             v_adj_snap=vy;
%             v_fw_snap = vy_fw;
%             interact = interaction.rho.y;
%             Kaa = K.rho.SH;
%             direction = 'Y';
%             prc=99.985;
%         elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%             v_adj_snap=vx;
%             v_fw_snap = vx_fw;
%             interact = interaction.rho.x;
%             Kaa = K.rho.x;
%             direction = 'X';
%             prc=99.97;
%         end
%     elseif(row==2)
%         v_adj_snap=vz;
%         v_fw_snap = vz_fw;
%         interact = interaction.rho.z;
%         Kaa = K.rho.z;
%         direction = 'Z';
%         prc=99.97;
%     elseif(row==3)
%         v_adj_snap=vy;
%         v_fw_snap = vy_fw;
%         interact = interaction.rho.y;
%         Kaa = K.rho.SH;
%         direction = 'Y';
%         prc=99.99;
%     end
    
    
    %% plot forward field --------------------------------------------------
    
    subplot(nplots,ncols,4*(i-1)+1);
    cla;
    
    %- plot field 
    pcolor(X,Z,plotting(i).v_fw');
    hold on
    axis image
    shading interp
    
    %- colour scale
    % scale = 0.8*max(max(abs(snapshot)));
    scale=prctile(abs(plotting(1).v_fw(:)),99.97);
    caxis([-1*scale scale]);
%     colormap(flipud(cm)); % only one color map per figure allowed... argh
    
    %- labels etc.
    xlabel('x [m]');
    ylabel('z [m]');
    title(['forward velocity [m/s] (',plotting(i).dir,')']);
    
    % timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
    for bips = 1:5;
        plot_adjoint_src_rec;
        text(0.05*Lx,0.92*Lz,timestamp) ;
        text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ...
            'HorizontalAlignment','right')
    end

    hold off
    
    
    %% plot adjoint field ---------------------------------------------
%     max_v = max(plotting(1).v(:));
    subplot(nplots,ncols,4*(i-1)+2);
    cla;
    hold on
    
    %- plot field
    pcolor(X,Z,plotting(i).v');
    plot_adjoint_src_rec;
    axis image
    shading interp
    
    %- colour scale
    % scale = 0.8*max(abs(v_current(:)));
    scale=prctile(abs(plotting(i).v(:)),99.97);
    caxis([-scale scale]);
%     colormap(flipud(cm)); % only one colormap per figure allowed...
    
    %- labels etc.
    xlabel('x [m]');
    ylabel('z [m]');
    title(['adjoint velocity [m/s] (',plotting(i).dir,')']);
    
    %- time stamp & other text
    % timestamp=['t [s] = ',num2str(nt*dt-(n-plot_every)*dt)];
    for bips = 1:5;
    text(0.05*Lx,0.92*Lz,timestamp) ;
    text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ...
        'HorizontalAlignment','right')
    end
    
    
    
    %% plot interaction ----------------------------------------------------
    
    subplot(nplots,ncols,4*(i-1)+3);
    cla;
    hold on
    
    %- plot field
    pcolor(X,Z,plotting(i).interact');
    plot_adjoint_src_rec;
    axis image
    shading interp
    
    %- colour scale
    scale = 0.8*max(max(abs(plotting(i).interact)));
    caxis([-scale scale]);
%     colormap(flipud(cm));
    
    %- labels etc.
    xlabel('x [m]');
    ylabel('z [m]');
    title(['forw. \cdot adj. velocity (',plotting(i).dir,')']);
    
    %- timestamp and other text
    for bips = 1:5;
    text(0.05*Lx,0.92*Lz,timestamp) ;
    text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ...
        'HorizontalAlignment','right')
    end
    
    %% plot kernel ---------------------------------------------------------
    
    subplot(nplots,ncols,4*(i-1)+4);
    cla;
    hold on
    
    %- plot field
    pcolor(X,Z,plotting(i).kernel');
    plot_adjoint_src_rec;
    
    %- colour scale
    cmax = prctile(plotting(i).kernel(:),99.97);
    caxis([-cmax cmax]);
%     colormap(cm_model);
    
    % colorbar;
    
    axis image
    shading interp
    xlabel('x [m]');
    ylabel('z [m]');
    title([plotting(i).knlname,' kernel']);
    
    % timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
    for bips = 1:5;
    text(0.05*Lx,0.92*Lz,timestamp) ;
    text(0.95*Lx,0.92*Lz,['max = \pm', num2str(cmax)], ...
        'HorizontalAlignment','right')
    end
    
    
    
end

%% record movie -------------------------------------------------------
    
    if strcmp(make_movie_adj,'yes')
    
        if exist('movie_index','var')
            movie_index=movie_index+1;
        else
            movie_index=1;
        end
        
        M(movie_index)=getframe(gcf);
        
    end
    
    
    if any(savetimes == temps)
        filename = [snapshotfile,'_adjoint_t',num2str(temps),'.png'];
        disp(['saving at time ',num2str(temps),' with filename ',filename, '!'])
        print(gcf,'-dpng','-r800',filename)
    end
           

hold off
% pause(0.01)






    % experimenting with the caxis (0) -- this was not satisfactory
    % caxis([-0.8*max(max(abs(K))) 0.8*max(max(abs(K)))]);
    % if (n>0.90*nt)
    %     %             disp 'BANG!!'
    %     %             factor=(nt-n+5)/nt;
    %     %             if(n>0.93*nt)
    %     factor=0.07;
    %     %             end
    %     %             factor
    %     caxis([-factor*max(max(abs(K))) factor*max(max(abs(K)))]);
    % end
    %------------------
    % experimenting with the caxis (1) -- blanking out the src location
    % this possibly didn't work because I was blanking out receivers!!!!
    % K_no_origsrc = K;
    % for i=1:nsorig
    %     iks = orig_x_id(i);
    %     zet = orig_z_id(i);
    %     K_no_origsrc(iks-5:iks+5 , zet-5:zet+5) = 0;
    % end
    % cmax = max(max(abs(K_no_origsrc)));
    %------------------
    % experimenting with the caxis (2) -- use standard deviation of K
    % works with 10*std for now, don't know about other configurations though
    % cmax = min(max(Kaa(:)),10*std(Kaa(:)));
    %------------------
    % experimenting with the caxis (3) -- location of maximum value
    % could still use it like in experiment (1) blanking out the max region
    % [value, k] = max(K(:));
    % [ie, jee] = ind2sub(size(K), k);
    % max_loc = [X(1,ie),Z(jee)];
    % p = plot(X(1,ie),Z(jee),'kx');
    % set(p,'Color','green','LineWidth',2);
    %------------------
    % experimenting with the caxis (3.1) -- around location of maximum value
    % [value, k] = max(kernel(:));
    % [ie, jee] = ind2sub(size(kernel), k);
    % K_no_origsrc = kernel;
    % blank=min(size(X))/15;
    % K_no_origsrc( ie-blank:ie+blank , jee-blank:jee+blank ) = 0;
    %------------------
    % experimenting with the caxis (4) -- percentile