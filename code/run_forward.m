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

path(path,'propagation/');
path(path,'../input/');

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
fig_mod = plot_model(model, param_plot);


%% source & receiver initialisation

% sources and receivers in forward simulations ----------------------------
if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))

    %- time axis ----------------------------------------------------------
    
    t=0:dt:dt*(nt-1);
    
    %- compute indices for source & receiver locations --------------------

    [src_x_id,src_z_id,ns] = compute_indices(src_x,src_z,Lx,Lz,dx,dz);
    [rec_x_id,rec_z_id,n_receivers] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
    
    %- (make and) plot source time function ---------------------------------

%     if isnan(stf_in);
%         switch stf_type
%             case {'delta_bp', 'heaviside_bp'}
%                 stf = make_source_time_function(t,stf_type,f_min,f_max);
%             case 'ricker'
%                 stf = make_source_time_function(t,stf_type,tauw_0, tauw, tee_0);
%         end
%     elseif isnumeric(stf_in);
%         stf = stf_in;
%     else
%         error('PANIC! stf type not recognised it seems');
%     end
%     fig_stf = plot_source_time_function(t,stf);
% 
%     % add the same source time function to all the sources.    
% %     prefactor = 1.875e5;
% %     prefactor = 1.875e5 / dx / dz;
%     prefactor = 1.0 / dx / dz;
%     for i=1:ns
%         stf_all{i}.x = prefactor* stf.*stf_PSV(1)./norm(stf_PSV);  % x direction
%         stf_all{i}.y = prefactor* stf;                           % y direction
%         stf_all{i}.z = prefactor* stf.*stf_PSV(2)./norm(stf_PSV);  % z direction
%     end
%     
%     stf = stf_all;
end

%% check values

%test_input_compatibility;



% %% initialise seismograms ------------------------------------------------- % Change this to save velocity seismograms only!
% % only for forward
% if(strcmp(wave_propagation_type,'SH'))
%     v_rec.y=zeros(n_receivers,nt);
% elseif(strcmp(wave_propagation_type,'PSV'))
%     v_rec.x    =zeros(n_receivers,nt);
%     v_rec.z    =zeros(n_receivers,nt);
% elseif(strcmp(wave_propagation_type,'both'))
%     v_rec.x    =zeros(n_receivers,nt);
%     v_rec.y    =zeros(n_receivers,nt);
%     v_rec.z    =zeros(n_receivers,nt);
% end



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
end
if strcmp(simulation_mode,'forward')
    if(strcmp(wave_propagation_type,'SH'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            save('../output/v_forward','vy_forward','-v7.3');
            save('../output/u_forward','uy_forward','-v7.3');
        end
    elseif(strcmp(wave_propagation_type,'PSV'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            save('../output/v_forward','vx_forward','vz_forward', '-v7.3');
            save('../output/u_forward','ux_forward','uz_forward', '-v7.3');
        end
    elseif(strcmp(wave_propagation_type,'both'))
        if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
            if (strcmp(save_u_fw,'yes') && strcmp(save_v_fw,'yes') )
                save('../output/v_forward','vx_forward','vy_forward','vz_forward', '-v7.3');
                save('../output/u_forward','ux_forward','uy_forward','uz_forward', '-v7.3');
            end

        end
    end
end

% disp 'saving seismograms...'
% savename = ['../output/',project_name,'.v_rec.mat'];
% save(savename, 'v_rec', '-v7.3');





%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    disp(['storing movie: ',movie_file]);
    % profile needs to be 'Uncompressed AVI' on my (ubuntu) computer. 
    % (Matlab says that win7+/mac10.x+ can do mpeg-4)
    writerObj=VideoWriter(movie_file,'Uncompressed AVI');  
    writerObj.FrameRate = 10;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end

end



function [useGivenModel, GivenModel, stf] = checkargs(arg)

narg = length(arg);

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