% calculate the step length

function [v_try,t,steplnArray,misfitArray] = calculate_step_length(teststep, ...
                                      niter, ...
                                      currentMisfit, kernel, ...
                                      v_obs)
%== 1. Preparation ===========================================================

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_doffer;

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
[mu,rho,lambda]=define_material_parameters(nx,nz,11);



%- determine the number of steps we'll try and divide teststep into that nr

if (niter == 1)
    nsteps = 11;
elseif (niter > 1)
    nsteps = 3;
else
    error('your inversion iteration seems to be <1');
end

steplnArray = 0: 2*teststep/(nsteps-1) : 2*teststep;



%- set up array in which the misfits will be stored

misfitArray = zeros(nsteps,1);
misfitArray(1) = currentMisfit.x + currentMisfit.y + currentMisfit.z;




%== 2. Filtering kernels =====================================================

%- smooth the kernels so that no singularities will be formed
%  Essentially a low-pass filter. Theory in http://www.dspguide.com/ch24/1.htm
%
%-- make smoothing filter (gaussian)
filt = fspecial('gaussian',[5 5],2);   %  (fspecial is normally part of the image processing toolbox but I adapted
                                       %  the (probably exactly equivalent) Octave code. fspecial found in tools/)
%-- smoothe all 'total' kernels (in rho mu lambda parametrisation)
% for sname = fieldnames(kernel)'

bips = figure;
for sname = {'rho' 'mu' 'lambda'}
    figure(bips);
%     clf;
%     disp(sname)
%     disp(kernel.(sname{1}).total(800))
    Ksmooth.(sname{1}) = conv2(kernel.(sname{1}).total, filt, 'same');
    % see note on conv2 as image filter by Hannes Ovrén on Stackoverflow:
    % http://stackoverflow.com/questions/1737960/mean-filter-for-smoothing-images-in-matlab
    
    
    % plot to check if ok:
    
    cmax = prctile(abs(kernel.(sname{1}).total(:)),99.97);
    if strcmp(sname{1},'rho')
        subplot(2,1,1)
        knl_title = ['original ',sname{1}];
        plot_kernel(X,Z,kernel.(sname{1}).total,knl_title,'fixed',cmax,[1 0]);
        colorbar;
        subplot(2,1,2)
        knl_title = ['smoothed ',sname{1}];
        plot_kernel(X,Z,Ksmooth.(sname{1}),knl_title,'fixed',cmax,[1 0]);
        colorbar;
    end
%     pause(1.0);
    
    
end



% STILL DO:
%- check if we don't end up with negative velocities (other params that
%  can't be negative?)


%== 3. Calculating updates and misfits =====================================================
disp 'starting forward calculation...'

%- START LOOP
for ntry = 2:nsteps
    disp(['... now doing iteration number ',num2str(ntry)]);
    
    steptry = steplnArray(ntry);
    
    %- for each step, update model using step length
    
    update_rho = steptry * Ksmooth.rho;
    update_mu  = steptry * Ksmooth.mu;
    update_lambda = steptry * Ksmooth.lambda;
    
    % minus update because the kernel is the gradient pointing in the
    % uphill direction, while we want to obtain the minimum which is in the
    % downhill direction.
    rho = rho - update_rho;
    mu  = mu - update_mu;
    lambda = lambda - update_lambda;
    
    
    %- for each step, run forward
    
    [v_try,t,~,~,~,~] = run_forward_update(rho,mu,lambda);
    
    
    %- for each step, calculate the misfit using make_adstf --> adapt mk adstf!
    [~, misfit] = make_all_adjoint_sources(v_try,v_obs,t,'waveform_difference');
    
    %- for each step, save the misfit to the misfit array
    misfitArray(ntry) = misfit.x + misfit.y + misfit.z;

end


%- plot the misfit versus the step

end