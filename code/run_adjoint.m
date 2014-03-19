%==========================================================================
% run adjoint simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

function [K] = run_adjoint(kerneltype)%(vx_forward, vy_forward, vz_forward)
disp 'Welcome to the ultimate adjoint experience!!'
disp 'initialising...'

%==========================================================================
%% set paths and read input
%==========================================================================

% path(path,'helper_programmes/');
path(path,'../input/');
path(path,'../output/');
path(path,'propagation/');
path(path,'../tools/');

% read input parameters from the input file
input_parameters;
nt=5*round(nt/5);   % to make sure nt is a multiple of 5 -- needed because
                    % run_forward saves the forward field every 5 timesteps
                    % only, and thus we need to compare.

% adaptations to input parameters:
% (needed to get run_adjoint to work like it should)
simulation_mode='adjoint';                                                 
orig_src_x = src_x;
orig_src_z = src_z;

% make figures appear right on the screen
set_figure_properties;

% load matfile containing the stored forward field. Takes a LONG time.
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


% read adjoint source time functions from file

for n=1:ns          % loop over sources
    for dir= 1:3    % loop over directions 1,2,3 = x,y,z
        fid=fopen(['../input/sources/adjoint/src_' num2str(n) '_' num2str(dir)],'r');
        stf(dir,n,1:nt)=fscanf(fid,'%g',nt);
        
        % plotting the source time functions
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



%==========================================================================
%% RUN WAVEFIELD PROPAGATION
%==========================================================================

fig_adjoint = figure;

% wavefield propagation is executed backwards in time for the adjoint case
run_wavefield_propagation;

