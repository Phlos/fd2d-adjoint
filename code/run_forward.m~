%==========================================================================
% run forward simulation
%
% input:
% -------
% to be supplied in the file ../input/input_parameters.m
%
% output:
%--------
% u: displacement seismograms
% t: time axis
% rec_x: x-coordinate of receiver positions
% rec_z: z-coordinate of receiver positions
%
%==========================================================================

function [ux,uy,uz,t,rec_x,rec_z]=run_forward

disp 'initialising...'
%==========================================================================
%% set paths and read input
%==========================================================================

path(path,'propagation/');
path(path,'../input/');

input_parameters;
nt=5*round(nt/5);


set_figure_properties;  % i.e. size of the figures & location on screen

%==========================================================================
% initialise simulation
%==========================================================================

%% material and domain ----------------------------------------------------

[mu,rho,lambda]=define_material_parameters(nx,nz,model_type);               

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

plot_model;


%% source & receiver initialisation

% sources and receivers in forward simulations ('forward', 'forward_green') -
if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))

    %- time axis ----------------------------------------------------------
    
    t=0:dt:dt*(nt-1);
    
    %- compute indices for source & receiver locations --------------------

    [src_x_id,src_z_id,ns] = compute_indices(src_x,src_z,Lx,Lz,dx,dz);
    [rec_x_id,rec_z_id,n_receivers] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
    
    %- make and plot source time function ---------------------------------

    make_source_time_function;
    plot_source_time_function;                                              % add direction here! (and later on where DS is calculated, an if for in which direction the stf is)
    
    for i=1:ns
        stf_all(i,:) = stf;
        stf = stf_all;
    end
    
end


%% initialise seismograms ------------------------------------------------- % Change this to save velocity seismograms only!
% only for forward
if(strcmp(wave_propagation_type,'SH'))
    uy=zeros(n_receivers,nt);
    v_rec_y=zeros(n_receivers,nt);
elseif(strcmp(wave_propagation_type,'PSV'))
    ux=zeros(n_receivers,nt);
    uz=zeros(n_receivers,nt);
    v_rec_x    =zeros(n_receivers,nt);
    v_rec_z    =zeros(n_receivers,nt);
elseif(strcmp(wave_propagation_type,'both'))
    ux=zeros(n_receivers,nt);
    uz=zeros(n_receivers,nt);
    uy=zeros(n_receivers,nt);
    v_rec_x    =zeros(n_receivers,nt);
    v_rec_y    =zeros(n_receivers,nt);
    v_rec_z    =zeros(n_receivers,nt);
end



%==========================================================================
%% RUN THE ACTUAL WAVEFIELD PROPAGATION
%==========================================================================
run_wavefield_propagation;




%==========================================================================
%% output 
%==========================================================================
disp 'saving output...'
%- store time-reversed forward field --------------------------------------

if strcmp(simulation_mode,'forward')
    if(strcmp(wave_propagation_type,'SH'))
        save('../output/v_forward_total','vy_forward','-v7.3');
    elseif(strcmp(wave_propagation_type,'PSV'))
        save('../output/vx_forward','vx_forward','vz_forward','-v7.3');
%         save('../output/vz_forward','vz_forward');
    elseif(strcmp(wave_propagation_type,'both'))
        save('../output/v_forward','vx_forward','vy_forward','vz_forward','-v7.3');
    end
end


%- displacement seismograms -----------------------------------------------

if(strcmp(wave_propagation_type,'SH'))
    uy=cumsum(v_rec_y,2)*dt;
    ux=uy*0;
    uz=uy*0;
elseif(strcmp(wave_propagation_type,'PSV'))
    ux=cumsum(v_rec_x,2)*dt;
    uz=cumsum(v_rec_z,2)*dt;
    uy=ux*0;
elseif(strcmp(wave_propagation_type,'both'))
    uy=cumsum(v_rec_y,2)*dt;
    ux=cumsum(v_rec_x,2)*dt;
    uz=cumsum(v_rec_z,2)*dt;
end

%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    disp 'storing movie...'
    % profile needs to be 'Uncompressed AVI' on my (ubuntu) computer. 
    % (Matlab says that win7+/mac10.x+ can do mpeg-4)
    writerObj=VideoWriter(movie_file,'Uncompressed AVI');   
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end