% calculate the step length

function [step, fig_lnsrch , steplnArray_total, misfitArray_total ] = calculate_step_length(teststep, niter, ...
                                      currentMisfit, misfit_init, ...
                                      Model_prev, K_abs, g_obs, sEventInfo, sEventObs, Model_start)
%% == 1. Preparation ======================================================

% % obtain useful parameters from input_parameters
% [~, ~, ~, use_grav, ~, ~, ~, ~, ~, ~] = get_input_info;

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

Model_bg = update_model(bg_model_type);

%- plot starting model
% fig_mod_prev = plot_model_diff(Model_prev, Model_bg, param_plot);


%- determine the number of steps we'll try and divide teststep by that nr
if (niter == 1)
    nsteps = 5;
elseif (niter > 1)
    nsteps = 3;
else
    error('your inversion iteration seems to be <1');
end

%- determine the step lengths to be tested: 0, teststep, 2* teststep
teststep = abs(teststep);
steplnArray = 0: 2*teststep/(nsteps-1) : 2*teststep;
disp(['number of step lengths we will investigate: ',num2str(nsteps), ...
      ' -- test step length ', num2str(teststep,'%3.1e')]);


%- set up array in which the misfits will be stored
misfitArray.total = zeros(1,nsteps);
misfitArray.total(1) = currentMisfit;



%% == 2. Calculating updates and misfits ==================================


%- START LOOP
for ntry = 2:nsteps
        
    steptry = steplnArray(ntry);
    disp ' ';
    disp '----------------------------------------------';
%     disp ' ';
    disp(['Step ',num2str(ntry), ' of ', num2str(nsteps), ...
          ' --- step length ', num2str(steptry,'%3.1e')]);

    misfitArray.total(ntry) = calc_misfit_perstep(K_abs, steptry, Model_prev, ...
        misfit_init, g_obs, sEventInfo, sEventObs, Model_start);

    %- give step nr, step length and misfit
    disp(['Step ',num2str(ntry), ...
          ': step length ', num2str(steptry,'%3.1e'), ...
          ' and misfit ', num2str(misfitArray.total(ntry)), ...
          ' (diff ', num2str((misfitArray.total(ntry) - currentMisfit) / currentMisfit ),')']);

end



% STILL DO:
%- check if we don't end up with negative velocities (other params that
%  can't be negative?)





%% == 3. Calculate the real step length ===================================

close all;
[step, minval, p, RSS, fig_lnsrch] = fit_polynomial(steplnArray, misfitArray.total);
FitGoodness = RSS / (max(misfitArray.total) - min(misfitArray.total));
% print(fig_lnsrch,'-dpng','-r400','../output/fig_linesearch.initial.png');
disp ' ';
disp(['Resulting step: step length ', num2str(step,'%3.1e'), ...
    ' --- misfit      ', num2str(minval,'%3.1e')]);


%% == 4. Catch bad polynomials ============================================
%- make a catch for:
%  max extremum 
%  negative step length
%  poorly fitting polynomial

% FitGoodness = S.normr / (max(misfitArray.total) - min(misfitArray.total));
% FitGoodness = RSS / (max(misfitArray.total) - min(misfitArray.total));

steplnArray_total = steplnArray;
misfitArray_total = misfitArray;

steplnArray_prev = steplnArray;
misfitArray_prev = misfitArray;

nextra = 1;

% make new test points closer or farther away from starting point depending
% on why the polyfit doesn't give a good step length
while ( (p(1) < 0 || step < 0 || FitGoodness > 0.1) && nextra <= 5 )
    disp ' ';
    disp '----------------------------------------------';
    disp 'we obtained a max extremum or the step is negative'
    if (p(1) < 0) && step < 0
        errorreason = 'neg_max';
        disp 'it was a maximum AND negative step size'
    elseif (p(1) < 0)
        errorreason = 'max';
        disp 'it was a maximum'
    elseif step < 0
        errorreason = 'neg';
        disp 'the step size was negative'
    elseif (FitGoodness > 0.1)
        errorreason = 'poorfit';
        disp(['the polyfit is poor: ', num2str(FitGoodness)]);
    else
        error('none of the above wtf???')
    end
    
    % determine new steplnarray and misfitarray
    if strcmp(errorreason,'neg_max')
%         max(steplnArray_prev(
        steplnArray_new = [0 , max(steplnArray_prev(:)), 2*max(abs(steplnArray_prev(:))) ];
        misfitArray_new.total(1) = misfitArray_prev.total(1);
        misfitArray_new.total(2) = misfitArray_prev.total(end);
        idxEmpty = 3;
%     elseif ( nextra ==1 && steplnArray_prev(2) > 2*stepInit)
%         steplnArray_new = [0 , stepInit , steplnArray_prev(2)];
%         misfitArray_new.total(1) = misfitArray_prev.total(1);
%         misfitArray_new.total(3) = misfitArray_prev.total(2);
%         idxEmpty = 2;
    else
        steplnArray_new = [0 , (1.0/3)*steplnArray_prev(2) , steplnArray_prev(2)];
        misfitArray_new.total(1) = misfitArray_prev.total(1);
        misfitArray_new.total(3) = misfitArray_prev.total(2);
        idxEmpty = 2;
    end
    teststep_new = steplnArray_new(idxEmpty);
    
    disp ' ';
    disp(['Testing NEW step ',num2str(nextra), ...
          ' --- step length ', num2str(teststep_new,'%3.1e')]);
    
    % calculate misfit for the new teststep
    [misfit_new, misfitseis, misfitgrav] = calc_misfit_perstep(K_abs, teststep_new, ...
        Model_prev, misfit_init, g_obs, sEventInfo, sEventObs, Model_start);
    disp ' ';
        %- give step nr, step length and misfit
    disp({['Extra step ',num2str(nextra), ...
          ': step length ', num2str(teststep_new,'%3.1e'), ...
          ' and misfit ', num2str(misfit_new,'%3.1e'), ...
          ' (diff) ', num2str((misfit_new - currentMisfit) / currentMisfit )];
          ['                  seis: ', num2str(misfitseis,'%3.1e')]; ...
          ['                  grav: ', num2str(misfitgrav,'%3.1e')]; ...
          });
      
      % save new misfit to appropriate place in new array and total array
    misfitArray_new.total(idxEmpty) = misfit_new;
    misfitArray_total.total = [misfitArray_total.total, misfit_new];
    steplnArray_total = [steplnArray_total, teststep_new];
    
    
    % determine whether we now get a minimum at positive step length
close all
    [step, minval, p, RSS, fig_lnsrch] = fit_polynomial(steplnArray_new, misfitArray_new.total);
    FitGoodness = RSS / (max(misfitArray.total) - min(misfitArray.total));
%     print(fig_lnsrch,'-dpng','-r400',['../output/fig_linesearch.extra-',num2str(nextra),'.png']);
    
    disp ' ';
    disp(['Found step: step length ', num2str(step,'%3.1e'), ...
        ' --- misfit ', num2str(minval,'%3.1e')...
        ' --- (diff ', num2str((minval-currentMisfit)/currentMisfit,'%3.1e'),')']);

    % if so, update step and p(1) (this will exit the loop)
    
    % add one to nextra
    nextra = nextra +1;
    
    % make the new 'prev' arrays
    steplnArray_prev = steplnArray_new;
    misfitArray_prev = misfitArray_new;
    

end

%% == 5. Output ===========================================================

%- do something with steplnArray_total and misfitArray_total here!
close all;
[~, ~, ~, ~, fig_lnsrch] = fit_polynomial(steplnArray_total, misfitArray_total.total);
fig_lnsrch.Children.Children(1).Color = 'w';
fig_lnsrch.Children.Children(2).LineStyle = '--';
% add the final polyfit and step to it
iksmin = min(steplnArray_total); iksmaks = max(steplnArray_total); 
iks = iksmin:(iksmaks-iksmin)/30:iksmaks;
ei = p(1)*iks.^2 + p(2)*iks + p(3);
hold on; echtefit = plot(iks,ei,'-r', 'LineWidth', 1); %echtefit.LineWidth = 2;
% calculate extremum value and plot
echtestap = plot(step,minval,'kx', 'LineWidth', 2); %echtestap.LineWidth = 2;



%% == 6. If linesearch still sucks, quit ==================================
    
%- check if extremum is a minimum, else error and quit
if (p(1) < 0) && step < 0
        disp 'Your step length line search obtained a NEGATIVE maximum!!'
elseif (p(1) < 0)
    error('Your step length line search obtained a maximum...');
elseif step < 0
    warning('The obtained step length is negative!!');
end
    
disp ' ';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
% disp '======================================';
disp(['Step length found: ',num2str(step,'%3.2e'), ...
      ' with misfit of ~', num2str(minval,'%3.2e'), ...
      ' (diff ', num2str((minval - currentMisfit) / currentMisfit ),')']);
% disp '======================================';
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
disp '----------------------------------------------------';
disp ' ';


end


%% subfunctions

function [misfittotal, misfitseis, misfitgrav] = ...
                       calc_misfit_perstep(K_abs, steptry, Model_prev, ...
                       misfit_init, g_obs, sEventInfo, sEventObs, Model_start)

% paths etc.

path(path,'../code');
path(path,'../input');
input_parameters;
set_figure_properties_bothmachines;

    %% calculate updated model using steptry
%     disp 'now updating model in calc_misfit_perstep'
    Model_prevfix = update_model(K_abs, steptry, Model_prev, parametrisation);
    if(strcmp(fix_velocities,'yes'))
        Model_try = fix_vs_vp(Model_prevfix, Model_start);
    else
        Model_try = Model_prevfix;
    end

    % no negative density
    min_rh = min(Model_try.rho(:));
    if min_rh < 0
        error('density inside the test model is negative somewhere!!')
    end

    
    %- some output
    max_mu = max(Model_try.mu(:) - Model_prev.mu(:));
    max_rh = max(Model_try.rho(:) - Model_prev.rho(:));
    max_la = max(Model_try.lambda(:) - Model_prev.lambda(:));
    min_mu = min(Model_try.mu(:) - Model_prev.mu(:));
    min_rh = min(Model_try.rho(:) - Model_prev.rho(:));
    min_la = min(Model_try.lambda(:) - Model_prev.lambda(:));
    
    disp(['Max diff w prev model -- mu: ',num2str(max_mu,'%3.2e'), '   rho: ', num2str(max_rh), ...
          '   lambda: ', num2str(max_la,'%3.2e')]);
    disp(['Min diff w prev model -- mu: ',num2str(min_mu,'%3.2e'), '   rho: ', num2str(min_rh), ...
          '   lambda: ', num2str(min_la,'%3.2e')]);
    
    fig_mod = plot_model(Model_try);

    % NEW as of 29-4-2015
    [misfit.total, misfit.seis, misfit.grav] = calc_misfits(Model_try, ...
        g_obs, misfit_init.grav, sEventInfo, sEventObs, misfit_init.seis, ...
        'noplot','notext');
    misfittotal = misfit.total;
    misfitseis = misfit.seis;
    misfitgrav = misfit.grav;
    
    close(fig_mod);
    
    
%         %% display
%     disp({['    step length ', num2str(steptry,'%3.1e'), ...
%           ' and misfit ', num2str(misfittotal,'%3.1e'), ...
%           ' (diff) ', num2str((misfit_new - currentMisfit) / currentMisfit )];
%           ['                  seis: ', num2str(misfitseis,'%3.1e')]; ...
%           ['                  grav: ', num2str(misfitgrav,'%3.1e')]; ...
%           });

end

function [step, minval, p, RSS, fig_lnsrch] = fit_polynomial(steplnArray, misfitArraytotal)

%- determine teststep & nsteps
maxstep = max(steplnArray(:));
nsteps = length(steplnArray);

%- fit quadratic
[p] = polyfit(steplnArray,misfitArraytotal,2);
%- determine minimum of quadratic --> step length
step = -p(2)/(2*p(1));

% goodness of fit
iks = steplnArray;
ei  = p(1)*iks.^2 + p(2)*iks + p(3);
RSS = sum((misfitArraytotal - ei).^2);

%- plot the misfit versus the step
fig_lnsrch = figure;
hold on;

% plot misfits
plot(steplnArray,misfitArraytotal,'-bo');

% plot fitting polynomial
minstep = min([steplnArray, 0, step]);
maxstep = max(maxstep,step);
iks = 0: maxstep/(20*(nsteps-1)) : maxstep;
ei  = p(1)*iks.^2 + p(2)*iks + p(3);
plot(iks,ei,'r');

% calculate extremum value and plot
minval = p(1)*step.^2 + p(2)*step + p(3);
plot(step,minval,'kx');
hold off;

end