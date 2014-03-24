%==========================================================================
% run adjoint simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

function [K] = run_adjoint(stf,kerneltype) %(vx_forward, vy_forward, vz_forward)
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

% stf=zeros(3,ns,nt);
fig_adjoint_stf = figure;

% direction = zeros(1,2);

% read adjoint source time functions from file

% if (strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
%     if strcmp(adjoint_source_component,'x')
%         direction(1)=1;
%     elseif strcmp(adjoint_source_component,'z')
%         direction(1)=3;
%     else
%         error('the direction of the source component that you gave is invalid')
%     end
% end

% if (strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
%     direction(2)=2;
% end

figure(fig_adjoint_stf);
for n=1:ns          % loop over sources
    for dir= 1:3    % loop over directions 1,2,3 = x,y,z
%     for dir = direction;
%         disp(['reading direction ',num2str(dir),'.'])
%         fid=fopen(['../input/sources/adjoint/src_' num2str(n) '_' num2str(dir)],'r');
%         stf(dir,n,1:nt)=fscanf(fid,'%g',nt);
        
        % plotting the source time functions
        
        
        subplot(3,1,dir);
        dir
        thee=0:dt:nt*dt-dt;
        length(thee);
        oempa=reshape(stf(dir,n,:),1,nt);
        size(oempa);
        size(thee);
%         clf;
        plot(thee,oempa);
    end
end



%==========================================================================
%% RUN WAVEFIELD PROPAGATION
%==========================================================================

fig_adjoint = figure;

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

% wavefield propagation is executed backwards in time for the adjoint case
run_wavefield_propagation;





%==========================================================================
%% store output
%==========================================================================

disp 'storing kernels...'

%  if(strcmp(wave_propagation_type,'SH'))
        save('../output/kernels','K','-v7.3');
%         save('../output/v_rec', 'v_rec_y', '-v7.3');
%     elseif(strcmp(wave_propagation_type,'PSV'))
%         save('../output/v_forward','vx_forward','vz_forward', ...
%                                     'v_rec_x','v_rec_z','-v7.3');
%         save('../output/v_rec', 'v_rec_x','v_rec_z', '-v7.3');
%     elseif(strcmp(wave_propagation_type,'both'))
%         save('../output/v_forward','vx_forward','vy_forward','vz_forward',...
%                                    'v_rec_x','v_rec_y','v_rec_z','-v7.3');
%         save('../output/v_rec', 'v_rec_x','v_rec_y','v_rec_z', '-v7.3');
%     end

if strcmp(make_movie_adj,'yes')
    disp(['storing movie: ',movie_file_adj]);
    % profile needs to be 'Uncompressed AVI' on my (ubuntu) computer. 
    % (Matlab says that win7+/mac10.x+ can do mpeg-4)
    writerObj=VideoWriter(movie_file_adj,'Uncompressed AVI');   
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end

