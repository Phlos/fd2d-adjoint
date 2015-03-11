% calculate the step length

function [step, fig_lnsrch ] = calculate_step_length(teststep, niter, ...
                                      currentMisfit, misfit_seis, misfit_grav, ...
                                      Model_prev, K_abs, v_obs, g_obs)
%== 1. Preparation ===========================================================

% % obtain useful parameters from input_parameters
% [~, ~, ~, use_grav, ~, ~, ~, ~, ~, ~] = get_input_info;

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

%- plot starting model
fig_mod_prev = plot_model(Model_prev, parametrisation);


%- determine the number of steps we'll try and divide teststep by that nr
if (niter == 1)
    nsteps = 3;
elseif (niter > 1)
    nsteps = 3;
else
    error('your inversion iteration seems to be <1');
end

%- determine the step lengths to be tested: 0, teststep, 2* teststep
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
    disp '----------------------------------------------';
%     disp ' ';
    disp(['Step ',num2str(ntry), ' of ', num2str(nsteps), ...
          ' --- step length ', num2str(steptry,'%3.1e')]);

    misfitArray.total(ntry) = calc_misfit_perstep(K_abs, steptry, Model_prev, ...
        misfit_grav, misfit_seis, g_obs, v_obs);

    %- give step nr, step length and misfit
    disp(['Step ',num2str(ntry), ...
          ': step length ', num2str(steptry,'%3.1e'), ...
          ' and misfit ', num2str(misfitArray.total(ntry))]);

end



% STILL DO:
%- check if we don't end up with negative velocities (other params that
%  can't be negative?)





%== 4. Calculate the real step length =====================================

[step, minval, p, fig_lnsrch] = fit_polynomial(steplnArray, misfitArray.total);


%== 5. Catch a maximum in the polynomial or a negative step length ========



%- make a catch for max extremum and negative step length

% teststep_prev = teststep;
steplnArray_prev = steplnArray;
misfitArray_prev = misfitArray;
nextra = 1;

% while the poly gives a max OR the steplength is negative
% keep using the smallest previous test value to see if we end up with a
% minimum / a positive step length etc.
while ( (p(1) < 0 || step < 0) && nextra < 5 )
    disp ' ';
    disp '----------------------------------------------';
    disp 'we obtained a max extremum or the step is negative'
    if (p(1) < 0)
        disp 'it was a maximum'
    elseif step < 0
        disp 'the step size was negative'
    else
        error('none of the above wtf???')
    end
    
    % determine new steplnarray
    if nextra ==1
        steplnArray_new = [0 , stepInit , steplnArray_prev(2)];
    else
        steplnArray_new = [0 , 0.5*steplnArray_prev(2) , steplnArray_prev(2)];
    end
    teststep_new = steplnArray_new(2);
    
    disp ' ';
    disp(['Testing NEW step ',num2str(nextra), ...
          ' --- step length ', num2str(teststep_new,'%3.1e')]);
    
    % determine new misfitarray
    misfitArray_new.total(1) = misfitArray_prev.total(1);
    misfitArray_new.total(3) = misfitArray_prev.total(2);
    
    % fill in value 2 of new misfit array
    misfitArray_new.total(2) = calc_misfit_perstep(K_abs, teststep_new, ...
        Model_prev, misfit_grav, misfit_seis, g_obs, v_obs);
    disp ' ';
    disp(['Result:          at step length ', num2str(teststep_new,'%3.1e'), ...
        ' misfit : ', num2str(misfitArray_new.total(2),'%3.1e')]);
    
    % determine whether we now get a minimum at positive step length
%     disp 'fitting the NEW polynomial'
    close(fig_lnsrch);
    [step, minval, p, fig_lnsrch] = fit_polynomial(steplnArray_new, misfitArray_new.total);
    
    
    disp ' ';
    disp(['Found step:    ',num2str(nextra), ...
        ' --- step length ', num2str(step,'%3.1e'), ...
        ' --- misfit      ', num2str(minval,'%3.1e')]);

    % if so, update step and p(1) (this will exit the loop)
    
    % add one to nextra
    nextra = nextra +1;
    
    % make the new 'prev' arrays
%     disp 'making the NEW ''prev'' arrays'
    steplnArray_prev = steplnArray_new;
    misfitArray_prev = misfitArray_new;
    
    % make something that stores all the extra misfits & tried steps into
    % one array that can be plotted after we're all ready and done. (needs
    % to be sorted as well)
end

%- check if extremum is a minimum, else error and quit
if (p(1) < 0)
    error('Your step length line search obtained a maximum...');
elseif step < 0
    warning('The obtained step length is negative!!');
end
    
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
% disp '======================================';
disp(['Step length found: ',num2str(step,'%3.2e'), ...
      ' with misfit of ~', num2str(minval,'%3.2e')]);
% disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '----------------------------------------------------';
disp ' ';




figure(fig_mod_prev);
pause(10);

close(fig_mod_prev);
close all;

end

function misfittotal = calc_misfit_perstep(K_abs, steptry, Model_prev, misfit_grav, misfit_seis, g_obs, v_obs)

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

    %% calculate updated model using steptry
    disp 'now updating model in calc_misfit_perstep'
    Model_try = update_model(K_abs, steptry, Model_prev, parametrisation);
    
    %- some output
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
      
    %% gravity misfit
    
    if strcmp(use_grav,'yes')
      [g_try, fig_grav] = calculate_gravity_field(Model_try.rho, rec_g);
      close(fig_grav);
      
      %- calculate gravity misfit:
      scaling_g = misfit_grav(1).total;
      [~, misf_g_test] = make_gravity_sources(g_try, g_obs, scaling_g);
      
      misfit.grav = misf_g_test.normd;
    end
    
    %% seismic misfit
    
    %- for each step, run forward update
    [v_try,t,~,~,~,~] = run_forward(Model_try);
%     close(gcf);
%     close(gcf);
%     close(gcf);
    
    %- for each step, calculate the misfit
%     if strcmp(normalise_misfits, 'byfirstmisfit')
%         scaling_s = misfit_seis{1}.total;
%     else
%         scaling_s = 1;
%     end
    [~, misf_s_test] = calc_misfitseis_adstf('waveform_difference',t,v_try,v_obs);
%     [~, misf_s_test] = make_all_adjoint_sources(v_try,v_obs,t, ...
%         'waveform_difference','auto', scaling_s);
    if ( strcmp(normalise_misfits, 'byfirstmisfit') )
        misf_s_test.normd = misf_s_test.total / ...
            misfit_seis{1}.total;
    else
        misf_s_test.normd = misf_s_test.total;
    end
    misfit.seis = misf_s_test.normd;


    %% combine the misfits
    if strcmp(use_grav,'yes')
        misfit.total = misfit.seis + misfit.grav;
    else
        misfit.total =  misfit.seis;
    end
    misfittotal = misfit.total;
    
    if strcmp(use_grav,'yes')
        disp(['seismic misfit:  ', ...
            num2str(misfit.seis,'%3.2e')])
        disp(['gravity misfit:  ', ...
            num2str(misfit.grav,'%3.2e')])
    end

end

function [step, minval, p, fig_lnsrch] = fit_polynomial(steplnArray, misfitArraytotal)

%- determine teststep & nsteps
maxstep = max(steplnArray(:));
nsteps = length(steplnArray);

%- fit quadratic
p = polyfit(steplnArray,misfitArraytotal,2);
%- determine minimum of quadratic --> step length
step = -p(2)/(2*p(1));

%- plot the misfit versus the step
fig_lnsrch = figure;
hold on;

% plot misfits
plot(steplnArray,misfitArraytotal,'b');

% plot fitting polynomial
iks = 0: maxstep/(20*(nsteps-1)) : maxstep;
ei  = p(1)*iks.^2 + p(2)*iks + p(3);
plot(iks,ei,'r');

% calculate extremum value and plot
minval = p(1)*step.^2 + p(2)*step + p(3);
plot(step,minval,'kx');
hold off;

end