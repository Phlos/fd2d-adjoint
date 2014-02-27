% plot the adjoint fields
%==========================================================================

figure(fig_adjoint);
ncols=4;
K=0;

if(strcmp(wave_propagation_type,'SH'))
    set(fig_adjoint,'OuterPosition',pos_adj_1)
    nrows=1;
elseif(strcmp(wave_propagation_type,'PSV'))
    set(fig_adjoint,'OuterPosition',pos_adj_2)
    nrows=2;
elseif(strcmp(wave_propagation_type,'both'))
    set(fig_adjoint,'OuterPosition',pos_adj_3)
    nrows=3;
end


%% every row is a direction: x, y or z
for row=1:nrows
    if(row==1)
        if(strcmp(wave_propagation_type,'SH'))
            v_current=vy;
            snapshot = vy_forward_snapshot;
            interaction = interaction_y;
            K = Ky;
            direction = 'Y';
        elseif(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
            v_current=vx;
            snapshot = vx_forward_snapshot;
            interaction = interaction_x;
            K = Kx;
            direction = 'X';
        end
    elseif(row==2)
        v_current=vz;
        snapshot = vz_forward_snapshot;
        interaction = interaction_z;
        K = Kz;
        direction = 'Z';
    elseif(row==3)
        v_current=vy;
        snapshot = vy_forward_snapshot;
        interaction = interaction_y;
        K = Ky;
        direction = 'Y';
    end

%% plot adjoint field ---------------------------------------------
max_v = max(v_current(:));
subplot(nrows,ncols,4*(row-1)+1);
cla;
pcolor(X,Z,v_current');
hold on
plot_adjoint_src_rec;

caxis([-0.8*max(max(abs(v_current))) 0.8*max(max(abs(v_current)))]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['adjoint velocity field [m/s] (',direction,' component)']);

timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;


%% plot forward field --------------------------------------------------

subplot(nrows,ncols,4*(row-1)+2);
cla;
pcolor(X,Z,snapshot');
hold on
plot_adjoint_src_rec;

caxis([-0.8*max(max(abs(snapshot))) 0.8*max(max(abs(snapshot)))]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['forward velocity field [m/s] (',direction,' component)']);

timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;

%% plot interaction ----------------------------------------------------

subplot(nrows,ncols,4*(row-1)+3);
cla;
hold on
pcolor(X,Z,interaction');
plot_adjoint_src_rec;

caxis([-0.8*max(max(abs(interaction))) 0.8*max(max(abs(interaction)))]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['interaction (forward \cdot adjoint) (',direction,' component)']);

timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;

%% plot kernel ---------------------------------------------------------

subplot(nrows,ncols,4*(row-1)+4);
cla;
hold on
pcolor(X,Z,K');
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
%     for j=iks-5:iks+5
%         for k=zet-5:zet+5
%             K_no_origsrc(j,k) = 0;
%         end
%     end
% end
% cmax = max(max(abs(K_no_origsrc)));
%------------------
% experimenting with the caxis (2) -- use standard deviation of K
% works with 10*std for now, don't know about other configurations though
cmax = min(max(K(:)),10*std(K(:)));
%------------------
% experimenting with the caxis (3) -- location of maximum value
% could still use it like in experiment (1) blanking out the max region
% [value, k] = max(K(:));
% [ie, jee] = ind2sub(size(K), k);
% max_loc = [X(1,ie),Z(jee)];
% p = plot(X(1,ie),Z(jee),'kx');
% set(p,'Color','green','LineWidth',2);
%------------------



caxis([-cmax cmax]);
colormap(cm);
axis image
shading interp
xlabel('x [m]');
ylabel('z [m]');
title(['sensitivity kernel (',direction,' component)']);

timestamp=['t [s] = ',num2str(nt*dt-(n-5)*dt)];
text(0.05*Lx,0.92*Lz,timestamp) ;


end
hold off
% pause(0.01)

