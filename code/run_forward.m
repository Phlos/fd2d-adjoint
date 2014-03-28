%==========================================================================
% run forward simulation
%
% input:
% -------
% to be supplied in the file ../input/input_parameters.m
%
% output:
%--------
% vx,vy,vz: velocity seismograms
% t: time axis
% rec_x: x-coordinate of receiver positions
% rec_z: z-coordinate of receiver positions
%
%==========================================================================

function [v_rec,t,u_fw,v_fw,rec_x,rec_z]=run_forward

disp 'initialising...'
%==========================================================================
%% set paths and read input
%==========================================================================

path(path,'propagation/');
path(path,'../input/');

input_parameters;
% wave_propagation_type
nt=5*round(nt/5);


set_figure_properties;  % i.e. size of the figures & location on screen

%==========================================================================
%% initialise simulation
%==========================================================================

%% material and domain ----------------------------------------------------

[mu,rho,lambda]=define_material_parameters(nx,nz,model_type);               

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

plot_model;


%% source & receiver initialisation

% sources and receivers in forward simulations ----------------------------
if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))

    %- time axis ----------------------------------------------------------
    
    t=0:dt:dt*(nt-1);
    
    %- compute indices for source & receiver locations --------------------

    [src_x_id,src_z_id,ns] = compute_indices(src_x,src_z,Lx,Lz,dx,dz);
    [rec_x_id,rec_z_id,n_receivers] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
    
    %- make and plot source time function ---------------------------------

    make_source_time_function;
    plot_source_time_function;
    
    % add the same source time function to all the sources.
    stf_all=zeros(3,ns,nt);
    
    for i=1:ns
        stf_all(1,i,:) = stf.*stf_PSV(1)./norm(stf_PSV);  % x direction
        stf_all(2,i,:) = stf;                           % y direction
        stf_all(3,i,:) = stf.*stf_PSV(2)./norm(stf_PSV);  % z direction
        stf = stf_all;
    end
    
end

%% check values

test_input_compatibility;

disp 'You need to clear big variables from your workspace - else i''m SLOW.'
clearok = input('''yes'' to exit script, ''no'' to continue. ','s');
if (strcmp(clearok,'no') || strcmp(clearok,'n'))
%     clearvars u_fw v_fw;
else
    error('Stopping script - do something with the variables yourself eh!')
end

%% initialise seismograms ------------------------------------------------- % Change this to save velocity seismograms only!
% only for forward
if(strcmp(wave_propagation_type,'SH'))
%     uy=zeros(n_receivers,nt);
    v_rec.y=zeros(n_receivers,nt);
elseif(strcmp(wave_propagation_type,'PSV'))
%     ux=zeros(n_receivers,nt);
%     uz=zeros(n_receivers,nt);
    v_rec.x    =zeros(n_receivers,nt);
    v_rec.z    =zeros(n_receivers,nt);
elseif(strcmp(wave_propagation_type,'both'))
%     ux=zeros(n_receivers,nt);
%     uz=zeros(n_receivers,nt);
%     uy=zeros(n_receivers,nt);
    v_rec.x    =zeros(n_receivers,nt);
    v_rec.y    =zeros(n_receivers,nt);
    v_rec.z    =zeros(n_receivers,nt);
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

u_fw.x = ux_forward;
u_fw.y = uy_forward;
u_fw.z = uz_forward;

v_fw.x = vx_forward;
v_fw.y = vy_forward;
v_fw.z = vz_forward;

% % This is taking way too long. Only compute if explicitly asked

if strcmp(simulation_mode,'forward')
    if(strcmp(wave_propagation_type,'SH'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            disp 'v_forward...'
            save('../output/v_forward','vy_forward','-v7.3');
            disp 'u_forward...'
            save('../output/u_forward','uy_forward','-v7.3');
        end
    elseif(strcmp(wave_propagation_type,'PSV'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            disp 'v_forward...'
            save('../output/v_forward','vx_forward','vz_forward', '-v7.3');
            disp 'u_forward...'
            save('../output/u_forward','ux_forward','uz_forward', '-v7.3');
        end
    elseif(strcmp(wave_propagation_type,'both'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
                disp 'v_forward...'
                save('../output/v_forward','vx_forward','vy_forward','vz_forward', '-v7.3');
                disp 'u_forward...'
                save('../output/u_forward','ux_forward','uy_forward','uz_forward', '-v7.3');
            end

        end
    end
end

disp 'saving seismograms...'
save('../output/v_rec', 'v_rec', '-v7.3');

copyfile('../input/input_parameters.m','../output/input_parameters_out.m')



%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    disp(['storing movie: ',movie_file]);
    % profile needs to be 'Uncompressed AVI' on my (ubuntu) computer. 
    % (Matlab says that win7+/mac10.x+ can do mpeg-4)
    writerObj=VideoWriter(movie_file,'Uncompressed AVI');   
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end