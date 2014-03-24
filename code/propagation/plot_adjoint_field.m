% plot the adjoint fields
%==========================================================================

figure(fig_adjoint);
ncols=4;
% K=0;



timestamp=['t [s] = ',num2str(nt*dt-(n-plot_every)*dt)];

disp(timestamp);

%% every row is a direction: x, y or z
for row=1:nrows
    if(row==1)
        if(strcmp(wave_propagation_type,'SH'))
            v_current=vy;
            snapshot = vy_fw_snapshot;
            interaction = interaction_vy;
            Kaa = K.rho.SH;
            direction = 'Y';
            prc=99.985;
        elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
            v_current=vx;
            snapshot = vx_fw_snapshot;
            interaction = interaction_vx;
            Kaa = K.rho.x;
            direction = 'X';
            prc=99.97;
        end
    elseif(row==2)
        v_current=vz;
        snapshot = vz_fw_snapshot;
        interaction = interaction_vz;
        Kaa = K.rho.z;
        direction = 'Z';
        prc=99.97;
    elseif(row==3)
        v_current=vy;
        snapshot = vy_fw_snapshot;
        interaction = interaction_vy;
        Kaa = K.rho.SH;
        direction = 'Y';
        prc=99.99;
    end

%% plot adjoint field ---------------------------------------------
max_v = max(v_current(:));
subplot(nrows,ncols,4*(row-1)+1);
cla;
pcolor(X,Z,v_current');
hold on
plot_adjoint_src_rec;

scale = 0.8*max(abs(v_current(:)));

caxis([-scale scale]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['adjoint velocity field [m/s] (',direction,' component)']);

% timestamp=['t [s] = ',num2str(nt*dt-(n-plot_every)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;
text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ... 
                      'HorizontalAlignment','right')


%% plot forward field --------------------------------------------------

subplot(nrows,ncols,4*(row-1)+2);
cla;
pcolor(X,Z,snapshot');
hold on
plot_adjoint_src_rec;

scale = 0.8*max(max(abs(snapshot)));

caxis([-scale scale]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['forward velocity field [m/s] (',direction,' component)']);

% timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;
text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ... 
                      'HorizontalAlignment','right')

%% plot interaction ----------------------------------------------------

subplot(nrows,ncols,4*(row-1)+3);
cla;
hold on
pcolor(X,Z,interaction');
plot_adjoint_src_rec;

scale = 0.8*max(max(abs(interaction)));

caxis([-scale scale]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['forward \cdot adjoint velocity (',direction,' component)']);

% timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;
text(0.95*Lx,0.92*Lz,['max = \pm', num2str(scale,'%3.1e')], ... 
                      'HorizontalAlignment','right')

%% plot kernel ---------------------------------------------------------

subplot(nrows,ncols,4*(row-1)+4);
cla;
hold on
pcolor(X,Z,Kaa');
plot_adjoint_src_rec;

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

cmax = prctile(Kaa(:),prc);
caxis([-cmax cmax]);

colormap(cm);

% colorbar;

axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['\rho kernel (',direction,' component)']);

% timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;
text(0.95*Lx,0.92*Lz,['max = \pm', num2str(cmax)], ... 
                      'HorizontalAlignment','right')


                  
       
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
           

hold off
% pause(0.01)

