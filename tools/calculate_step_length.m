% calculate the step length

function [step] = calculate_step_length(teststep, niter, ...
                                      currentMisfit, Params_prev, K_rel, ...
                                      v_obs)
%== 1. Preparation ===========================================================

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_doffer;

%- plot starting model
fig_mod_prev = plot_model(Params_prev);


%- determine the number of steps we'll try and divide teststep by that nr
if (niter == 1)
    nsteps = 5;
elseif (niter > 1)
    nsteps = 3;
else
    error('your inversion iteration seems to be <1');
end

steplnArray = 0: 2*teststep/(nsteps-1) : 2*teststep;
disp(['number of step lengths we will investigate: ',num2str(nsteps), ...
      ' -- test step length ', num2str(teststep,'%3.1e')]);


%- set up array in which the misfits will be stored
misfitArray = zeros(1,nsteps);
misfitArray(1) = currentMisfit.total;



% In the end, I decided to filter the kernels after all. But this is done
% (slightly inefficiently, because multiple times) inside the update_model
% script instead of elsewhere. I could move the kernel filtering over to
% calculate_other_kernels, I suppose... --- Nienke Blom, 5 July 2014
% Ksmooth = smooth_kernels(kernel,11);



%== 2. Calculating updates and misfits =====================================================


%- START LOOP
for ntry = 2:nsteps
        
    steptry = steplnArray(ntry);
    disp ' ';
    disp '==============================================';
%     disp ' ';
    disp(['Now testing step ',num2str(ntry), ' of ', num2str(nsteps), ...
          ' with step length ', num2str(steptry,'%3.1e')]);


    %- calculate updated model using steptry
    Params_try = update_model(K_rel, steptry, Params_prev);
    
    max_mu = max(Params_try.mu(:));
    max_rh = max(Params_try.rho(:));
    max_la = max(Params_try.lambda(:));
    min_mu = min(Params_try.mu(:));
    min_rh = min(Params_try.rho(:));
    min_la = min(Params_try.lambda(:));
    
    disp(['Maxima -- mu: ',num2str(max_mu,'%3.2e'), '   rho: ', num2str(max_rh), ...
          '   lambda: ', num2str(max_la,'%3.2e')]);
    disp(['Minima -- mu: ',num2str(min_mu,'%3.2e'), '   rho: ', num2str(min_rh), ...
          '   lambda: ', num2str(min_la,'%3.2e')]);
    
    
    
    %- for each step, run forward update
    [v_try,t,~,~,~,~] = run_forward_update(Params_try.rho, ...
                                           Params_try.mu, ...
                                           Params_try.lambda);
    close(gcf);
    close(gcf);
    
    %- for each step, calculate the misfit using make_adstf --> adapt mk adstf!
    [~, misfit] = make_all_adjoint_sources(v_try,v_obs,t,'waveform_difference','auto');
    
    %- for each step, save the misfit to the misfit array
    if strcmp(wave_propagation_type,'both')
        misfitArray(ntry) = misfit.x + misfit.y + misfit.z;
    elseif strcmp(wave_propagation_type,'PSV')
        misfitArray(ntry) = misfit.x + misfit.z;
    elseif strcmp(wave_propagation_type,'SH')
        misfitArray(ntry) = misfit.y;
    end
    

    disp(['Step ',num2str(ntry), ...
          ': step length ', num2str(steptry,'%3.1e'), ...
          ' and misfit ', num2str(misfitArray(ntry))]);


end



% STILL DO:
%- check if we don't end up with negative velocities (other params that
%  can't be negative?)





%== 4. Calculate the real step length =====================================

%- fit quadratic
p = polyfit(steplnArray,misfitArray,2);
%- determine minimum of quadratic --> step length
step = -p(2)/(2*p(1));

%- plot the misfit versus the step
fig_linesearch = figure;
hold on;

% plot misfits
plot(steplnArray,misfitArray,'b');

% plot fitting polynomial
iks = 0: teststep/(10*(nsteps-1)) : 2*teststep;
ei  = p(1)*iks.^2 + p(2)*iks + p(3);
plot(iks,ei,'r');

% calculate extremum value and plot
minval = p(1)*step.^2 + p(2)*step + p(3);
plot(step,minval,'kx');
hold off;

% check if extremum is a minimum, else error and quit
if (p(1) < 0)
    error('Your step length line search obtained a maximum...')
end

disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
% disp '======================================';
disp(['Step length found: ',num2str(step,'%3.2e'), ...
      ' with misfit of ~', num2str(minval,'%3.2e')]);
% disp '======================================';
disp ' ';

% this is not necessary: the negative gradient direction is ensured in
% update_model already
% step = -step;


% save figure
figname = '../output/iter.step-linesearch.png';
print(fig_linesearch,'-dpng','-r400',figname);
close(fig_linesearch);

close(fig_mod_prev);


end