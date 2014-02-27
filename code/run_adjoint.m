%==========================================================================
% run adjoint simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

function [Kx,Ky,Kz] = run_adjoint
disp 'Welcome to the ultimate adjoint experience!!'
disp 'initialising...'
%==========================================================================
%% set paths and read input
%==========================================================================
%         global dummy
% path(path,'helper_programmes/');
path(path,'../input/');
path(path,'../output/');
path(path,'propagation/');

% input parameters
input_parameters;       % read all the input parameters from the input file
nt=5*round(nt/5);

% adaptations to input parameters:
% (needed to get run_adjoint to work like it should)
simulation_mode='adjoint';                                                      % lichtelijk knullige manier om hem adjoint te laten werken

orig_src_x = src_x;
orig_src_z = src_z;


set_figure_properties;  % i.e. size of the figures & location on screen

disp '.... and loooooaaaading v_forward....'
load ../output/v_forward.mat


%==========================================================================
% initialise simulation
%==========================================================================

%% material and domain ----------------------------------------------------

[mu,rho,lambda]=define_material_parameters(nx,nz,model_type);
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

plot_model; % hoeft niet in adjoint? is tenslotte 't zelfde model


%% read adjoint sources ---------------------------------------------------
disp 'inserting adjoint sources...'
                    %     adjoint_sources = matfile('../input/sources/adjoint/adjoint_sources.mat');
                    %     if(not(exist('adjoint_stf')))
                    %         adjoint_stf=adjoint_sources.adjoint_stf; % could do this w/ load too
                    %     end
    
fid=fopen('../input/sources/adjoint/source_locations','r');
src_x=zeros(1);
src_z=zeros(1);

k=1;
while (feof(fid)==0)
    src_x(k)=fscanf(fid,'%g',1);
    src_z(k)=fscanf(fid,'%g',1);
    fgetl(fid);
    k=k+1;
end

fclose(fid);


%- compute indices for adjoint source locations -----------------------

[src_x_id,src_z_id,ns] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
[orig_x_id,orig_z_id,nsorig] = compute_indices(orig_src_x,orig_src_z,Lx,Lz,dx,dz);

%%
%- read adjoint source time functions + plot 'em --------------------------

stf=zeros(3,ns,nt);
fig_adjoint_stf = figure;

for n=1:ns          % loop over sources
    for dir= 1:3    % loop over directions 1,2,3 = x,y,z
        fid=fopen(['../input/sources/adjoint/src_' num2str(n) '_' num2str(dir)],'r');
        stf(dir,n,1:nt)=fscanf(fid,'%g',nt);
        
        % plotting
        thee=0:dt:nt*dt-dt;
        length(thee);
        oempa=reshape(stf(dir,n,:),1,nt);
        size(oempa);
        size(thee);
        figure(fig_adjoint_stf);
        clf;
        plot(thee,oempa);
        pause(0.1);
    end
end




%% initialise dynamic fields ----------------------------------------------
% both forward and adjoint
% initialise_dynamic_fields;  % this just makes all dynamic field (v, stress,
                            % derivatives of v and stress wherever needed
                            % with zeros(dimensions).


% initialise absorbing boundary taper a la Cerjan ------------------------

% absbound=ones(nx,nz);
% init_absbound;

%==========================================================================
%% iterate
%==========================================================================
% disp 'iterating....'
fig_adjoint = figure;

run_wavefield_propagation;
% time loop over the iterations. NOTE - in the adjoint case this goes
% BACKWARDS in time.
% for n=1:nt
% 
%     
%     %% compute divergence of current stress tensor and add external forces
%     % same as forward
%     
%     DS=div_s(sxy,szy,dx,dz,nx,nz,order);
% 
%     for i=1:ns
%         DS(src_x_id(i),src_z_id(i))=DS(src_x_id(i),src_z_id(i))+stf(i,n);
%     end
% 
%     
%     %% update velocity field ----------------------------------------------
%     % same as forward
%     vy=vy+dt*DS./rho;
%     vy=vy.*absbound;
%     
%     %- compute derivatives of current velocity and update stress tensor ---
%     
%     sxy=sxy+dt*mu(1:nx-1,1:nz).*dx_v(vy,dx,dz,nx,nz,order);
%     szy=szy+dt*mu(:,1:nz-1).*dz_v(vy,dx,dz,nx,nz,order);
%     
% 
%     
% end
