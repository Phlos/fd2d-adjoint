%==========================================================================
% run adjoint simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

function [K] = run_adjoint(u_fw,v_fw,stf,varargin)
% disp 'Welcome to the ultimate adjoint experience!!'
disp 'initialising...'

% check varargin (i.e. if parameters are to be taken from input_parameters
% or from a previous model run)
[updateParams, updateable] = checkargs(varargin(:));

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
set_figure_properties_doffer;

% extract the u_fw and v_fw fields
if(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    ux_forward = u_fw.x;
    uz_forward = u_fw.z;
    vx_forward = v_fw.x;
    vz_forward = v_fw.z;
end
if(strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both'))
    uy_forward = u_fw.y;
    vy_forward = v_fw.y;
end

%==========================================================================
% initialise simulation
%==========================================================================

%% material and domain ----------------------------------------------------

if strcmp(updateParams,'no')
    [mu,rho,lambda]=define_material_parameters(nx,nz,model_type);
elseif strcmp(updateParams,'yes')
    disp 'updating parameters...'
    rho = updateable.rho;
    mu = updateable.mu;
    lambda = updateable.lambda;
end

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

model.rho = rho;
model.mu = mu;
model.lambda = lambda;
fig_mod = plot_model(model, parametrisation); % hoeft niet in adjoint? is tenslotte 't zelfde model


%% read adjoint sources ---------------------------------------------------
% disp 'inserting adjoint sources...'
%     adjoint_sources = matfile('../input/sources/adjoint/adjoint_sources.mat');
%     if(not(exist('adjoint_stf')))
%         adjoint_stf=adjoint_sources.adjoint_stf; % could do this w/ load too
%     end

src_x = rec_x;
src_z = rec_z;


%- compute indices for adjoint source locations -----------------------

[src_x_id,src_z_id,ns] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
[orig_x_id,orig_z_id,nsorig] = compute_indices(orig_src_x,orig_src_z,Lx,Lz,dx,dz);

%%


%- read adjoint source time functions + plot 'em --------------------------

% fig_adjoint_stf = figure;

% figure(fig_adjoint_stf);
% for n=1:ns          % loop over sources
%     for dir= 1:3    % loop over directions 1,2,3 = x,y,z
% %         disp(['reading direction ',num2str(dir),'.'])
% %         fid=fopen(['../input/sources/adjoint/src_' num2str(n) '_' num2str(dir)],'r');
% %         stf(dir,n,1:nt)=fscanf(fid,'%g',nt);
%         
%         % plotting the source time functions
%         subplot(3,1,dir);
%         thee=0:dt:nt*dt-dt;
%         oempa=reshape(stf(dir,n,:),1,nt);
%         plot(thee,oempa);
%     end
%     pause(0.05);
% end



%==========================================================================
%% RUN WAVEFIELD PROPAGATION
%==========================================================================


% wavefield propagation is executed backwards in time for the adjoint case
run_wavefield_propagation;

% %% normalise the kernels
% 
% % disp(['normalising the seismic kernels by normfac = ', num2str(normfac,'%3.2e')]);
% params = fieldnames(K);
% for m = 1:length(params)
% %     disp(['i ', num2str(m),' PARAM: ', params{m}])
%     comps = fieldnames(K.(params{m}));
%     for l = 1:length(comps);
% %         disp(['j ', num2str(l),'   comp: ', comps{l}])
%         K.(params{m}).(comps{l}) = normfac * K.(params{m}).(comps{l});
%     end
% end

%==========================================================================
%% store output
%==========================================================================

% disp 'storing kernels...'
% kernelsavename = ['../output/',project_name,'.kernels.mat'];
% save(kernelsavename,'K','-v7.3');


if strcmp(make_movie_adj,'yes')
    disp(['storing movie: ',movie_file_adj]);
    % profile needs to be 'Uncompressed AVI' on my (ubuntu) computer. 
    % (Matlab says that win7+/mac10.x+ can do mpeg-4)
    writerObj=VideoWriter(movie_file_adj,'Uncompressed AVI');   
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end


% close(fig_adjoint_stf);
if(nt>plot_every)
    close(fig_adjoint);
end
close(fig_mod);

end


%% function check args
function [updateModel, Model] = checkargs(arg)

narg = length(arg);

if ( narg == 1 && isstruct(arg{1}) )
        
    % extract the updateable parameters
    updateModel = 'yes';
    Model = arg{1};

elseif narg == 0;
    updateModel = 'no';
    Model = NaN;
    disp(['number of input arguments to run_adjoint: ', num2str(narg)]);
    return
    
else
    error('ununderstood input to run_adjoint???')
    
end


end

% function [updateModel, Model, normfac] = checkargs(arg)
% 
% narg = length(arg);
% 
% if narg == 2
%     if (isstruct(arg{1}) && isnumeric(arg{2}))
%         updateModel = 'yes';
%         Model = arg{1};
%         normfac = arg{2};
%     else
%         error('input to run_adjoint not understood');
%     end
% 
% elseif narg == 1 
%     if isstruct(arg{1});
%         
%         % extract the updateable parameters
%         updateModel = 'yes';
%         Model = arg{1};
%         normfac = 1;
%     elseif isnumeric(arg{1})
%         updateModel = 'no';
%         Model = 0;
%         normfac = arg{1};
%     end
% 
% elseif narg == 0;
%     updateModel = 'no';
%     Model = 0;
%     normfac = 1;
%     disp(['number of input arguments to run_adjoint: ', num2str(narg)]);
%     return
%     
% else
%     error('more than two optional input to run_adjoint???')
%     
% end
% 
% 
% end
