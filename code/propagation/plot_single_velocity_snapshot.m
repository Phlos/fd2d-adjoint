%==========================================================================
% Function to take a snapshot of the velocity field at a given time.
%
% input:
% time     timestep at which the snapshot must be taken
%
% output:
% a plot with the velocity fields required. 
%
% NOTE: Care must be taken that the input_parameters correspond to the
% v_forward(s) that will be loaded from which the snapshot are taken! If
% they are not the same, the output is absolute and utter bogus. 
%==========================================================================

function plot_single_velocity_snapshot(time)
%%
% load the input_parameters and calculate the (rounded) 
% From this file, the size of the domain, the timestep and other parameters
% will be loaded.

path(path,'../code/');

input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

n=5*round((time/dt)/5); % this rounding needs to take place so that only 
                        % plots that are actually saved in the v_forward 
                        % are requested.
n_bigass=min(round((nt-n)/5+1),nt/5);
% time_rounded=n*dt;

vy=zeros(nx,nz);
vx=zeros(nx,nz);
vz=zeros(nx,nz);

set_figure_properties;
load cm_velocity;


%%

testing_v_forward = matfile('../output/v_forward.mat');
if(not(exist('vy_forward') && exist('vx_forward') && exist('vz_forward')))
    vx_snapshot=testing_v_forward.vx_forward(n_bigass,:,:);
    vy_snapshot=testing_v_forward.vy_forward(n_bigass,:,:);
    vz_snapshot=testing_v_forward.vz_forward(n_bigass,:,:);
    
    vx=reshape(vx_snapshot,nx,nz);
    vy=reshape(vy_snapshot,nx,nz);
    vz=reshape(vz_snapshot,nx,nz);
end

%% discarded bit that uses 'load' instead of 'matfile'
% if(not(exist('v_forward') && exist('vx_forward') && exist('vz_forward')))
%     disp('LADENNNN vy')
%     load ../output/v_forward.mat
%     disp('LADENNNN vx')
%     load ../output/vx_forward.mat
%     disp('LADENNNN vz')
%     load ../output/vz_forward.mat
% end

% % Calculate the location of the correct snapshot and load info in v, vx, vz
% % for n_bigass=1:round(nt/5)
% for i=1:nx
%     for j=1:nz
%         v(i,j)=v_forward(n_bigass,i,j);
%         vx(i,j)=vx_forward(n_bigass,i,j);
%         vz(i,j)=vz_forward(n_bigass,i,j);
%     end
% end


%%
% actual plotting
fig_vel = figure;                                                                                        %%%%%%%%%%% why not figure(number)? 
if(strcmp(wave_propagation_type,'both'))
set(fig_vel,'OuterPosition',pos_vel_nplots3);                                                                                                % then it plots only in that plot I think ---> nope
else
    set(fig_vel,'OuterPosition',pos_vel);  
end

plot_velocity_field;
% % end