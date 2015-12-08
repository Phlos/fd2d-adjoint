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

function [v_rec,t,u_fw,v_fw,rec_x,rec_z]=run_forward(varargin)

[updateParams, updateable, stf] = checkargs(varargin(:));

if exist('prevmsg') reverseStr = repmat(sprintf('\b'), 1, length(prevmsg));
else reverseStr = '';
end
prevmsg = sprintf('initialising...');
fprintf([reverseStr, prevmsg]); 
% disp 'initialising...'
%==========================================================================
%% set paths and read input
%==========================================================================

% path(path,'propagation/');
% path(path,'../input/');

input_parameters;

sfe = store_fw_every;
nt=sfe*round(nt/sfe);

set_figure_properties_bothmachines;  % i.e. size of the figures & location on screen



%==========================================================================
%% initialise simulation
%==========================================================================

%% material and domain ----------------------------------------------------

if strcmp(updateParams,'no')
    [mu,rho,lambda]=define_material_parameters(nx,nz,model_type);
elseif strcmp(updateParams,'yes')
    if exist('prevmsg') reverseStr = repmat(sprintf('\b'), 1, length(prevmsg));
    else reverseStr = '';
    end
    prevmsg = sprintf('updating parameters...');
    fprintf([reverseStr, prevmsg]);
%     disp 'updating parameters...'
    rho = updateable.rho;
    mu = updateable.mu;
    lambda = updateable.lambda;
end

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

model.rho = rho;
model.mu = mu;
model.lambda = lambda;
% fig_mod = plot_model(model, param_plot);


%% source & receiver initialisation

% sources and receivers in forward simulations ----------------------------
if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))

    %- time axis ----------------------------------------------------------
    
    t=0:dt:dt*(nt-1);
    
    %- compute indices for source & receiver locations --------------------

    [src_x_id,src_z_id,ns] = compute_indices(src_x,src_z,Lx,Lz,dx,dz);
    [rec_x_id,rec_z_id,n_receivers] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
    
end

%% check values

%test_input_compatibility;


%==========================================================================
%% RUN THE ACTUAL WAVEFIELD PROPAGATION
%==========================================================================
run_wavefield_propagation;

%==========================================================================
%% output 
%==========================================================================

%- store time-reversed forward field --------------------------------------

if(strcmp(wave_propagation_type,'PSV') || strcmp(wave_propagation_type,'both'))
    u_fw.x = ux_forward;
    u_fw.z = uz_forward;
    v_fw.x = vx_forward;
    v_fw.z = vz_forward;
end
if(strcmp(wave_propagation_type,'SH') || strcmp(wave_propagation_type,'both') )
    u_fw.y = uy_forward;
    v_fw.y = vy_forward;
end



% % This is taking way too long. Only compute if explicitly asked
if (strcmp(save_u_fw,'yes') || (strcmp(save_v_fw,'yes')) )
    disp 'saving u_fw, v_fw output to file...'
    save(['../output/forwardfield.mat'], 'u_fw', 'v_fw', '-v6');
end


%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    disp(['storing movie: ',movie_file]);
    % profile can now be 'Motion JPEG AVI' (before 'Uncompressed AVI')
    % much more efficient in output file size
    writerObj=VideoWriter(movie_file,'Motion JPEG AVI');  
    writerObj.FrameRate = 10;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
    % frame rate 20 fps:
    writerObj=VideoWriter([movie_file,'.20-fps'],'Motion JPEG AVI');
    writerObj.FrameRate = 20;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end

end



function [useGivenModel, GivenModel, stf] = checkargs(arg)

narg = length(arg);

input_parameters;


if narg == 2;
    useGivenModel = 'yes';
    GivenModel = arg{1};
    stf = arg{2};
elseif narg == 1;
    % extract the updateable parameters
    if isstruct(arg{1})
        useGivenModel = 'yes';
        GivenModel = arg{1};
        stf = NaN;
    elseif iscell(arg{1})
        stf = arg{1};
        useGivenModel = 'no';
        GivenModel = NaN;
    end
%     updatables = fieldnames(updateInput);
else
    useGivenModel = 'no';
    GivenModel = 0;
    stf = zeros(1,nt); stf(1) = 1; 
    warning('using dummy delta function stf!!!');
%     disp(['number of input arguments to run_forward: ', num2str(narg)]);
    return
    
end

% for i = 1:length(updatables)
%    disp(['now updating ',updatables{i} ]);
% %    disp([cell2str(updatables{i}), ' = ', cell2str(updatables{i}), ' + ', updateInput, '.', updatables{i}]);
%    updatables{i} = updatables{i} + updateInput.updatables{i};
%    disp(updatables{i});
% end


end