% calculate the step length

function [step, fig_linesearch ] = calculate_step_length(teststep, niter, ...
                                      currentMisfit, misfit_seis, misfit_grav, ...
                                      Model_prev, K_abs, v_obs, g_obs)
%== 1. Preparation ===========================================================

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

%- plot starting model
fig_mod_prev = plot_model(Model_prev);


%- determine the number of steps we'll try and divide teststep by that nr
if (niter == 1)
    nsteps = 3;
elseif (niter > 1)
    nsteps = 3;
else
    error('your inversion iteration seems to be <1');
end

steplnArray = 0: 2*teststep/(nsteps-1) : 2*teststep;
disp(['number of step lengths we will investigate: ',num2str(nsteps), ...
      ' -- test step length ', num2str(teststep,'%3.1e')]);


%- set up array in which the misfits will be stored
misfitArray.total = zeros(1,nsteps);
misfitArray.total(1) = currentMisfit;



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
    Model_try = update_model(K_abs, steptry, Model_prev, parametrisation);
    
    max_mu = max(Model_try.mu(:));
    max_rh = max(Model_try.rho(:));
    max_la = max(Model_try.lambda(:));
    min_mu = min(Model_try.mu(:));
    min_rh = min(Model_try.rho(:));
    min_la = min(Model_try.lambda(:));
    
    disp(['Maxima -- mu: ',num2str(max_mu,'%3.2e'), '   rho: ', num2str(max_rh), ...
          '   lambda: ', num2str(max_la,'%3.2e')]);
    disp(['Minima -- mu: ',num2str(min_mu,'%3.2e'), '   rho: ', num2str(min_rh), ...
          '   lambda: ', num2str(min_la,'%3.2e')]);
%     dada
      
    %% gravity
      
      [g_try, fig_grav] = calculate_gravity_field(Model_try.rho, rec_g);
      close(fig_grav);
      
      %- calculate gravity misfit:
      scaling_g = misfit_grav(1).total;
      [~, misf_g_test] = make_gravity_sources(g_try, g_obs, scaling_g);
      
      misfitArray.grav(ntry) = misf_g_test.normd;
      
    %% seismic
    
    %- for each step, run forward update
    [v_try,t,~,~,~,~] = run_forward(Model_try);
    close(gcf);
    close(gcf);
    close(gcf);
    
    %- for each step, calculate the misfit using make_adstf --> adapt mk adstf!
    scaling_s = misfit_seis(1).total;
    [~, misf_s_test] = make_all_adjoint_sources(v_try,v_obs,t, ...
        'waveform_difference','auto', scaling_s);
    misfitArray.seis(ntry) = misf_s_test.normd;
    
    % don't need to do this as misfit_s already has total misfit!!
%     %- for each step, save the misfit to the misfit array
%     if strcmp(wave_propagation_type,'both')
%         misfitArray.seis(ntry) = misfit_s.x + misfit_s.y + misfit_s.z;
%     elseif strcmp(wave_propagation_type,'PSV')
%         misfitArray.seis(ntry) = misfit_s.x + misfit_s.z;
%     elseif strcmp(wave_propagation_type,'SH')
%         misfitArray.seis(ntry) = misfit_s.y;
%     end


    misfitArray.total(ntry) = misfitArray.seis(ntry) + misfitArray.grav(ntry);

    
    disp(['gravity misfit:  ', ...
        num2str(misfitArray.seis(ntry),'%3.2e')])
    disp(['seismic misfit:  ', ...
        num2str(misfitArray.grav(ntry),'%3.2e')])
    disp(['Step ',num2str(ntry), ...
          ': step length ', num2str(steptry,'%3.1e'), ...
          ' and misfit ', num2str(misfitArray.total(ntry))]);


end



% STILL DO:
%- check if we don't end up with negative velocities (other params that
%  can't be negative?)





%== 4. Calculate the real step length =====================================

%- fit quadratic
p = polyfit(steplnArray,misfitArray.total,2);
%- determine minimum of quadratic --> step length
step = -p(2)/(2*p(1));

%- plot the misfit versus the step
fig_linesearch = figure;
hold on;

% plot misfits
plot(steplnArray,misfitArray.total,'b');

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






close(fig_mod_prev);


end