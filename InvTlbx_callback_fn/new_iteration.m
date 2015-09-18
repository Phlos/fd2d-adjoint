function [usr_par] = new_iteration( it, m, ModRandString, jm, gm, usr_par) 
% NEW_ITERATION This auxiliary function is called before a new iteration is started.
% You might use it for writing output files, plotting, ...
%
%
% Input: 
% m : current model
% it : current iteration number
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% usr_par : auxiliary user defined parameters (optional)



%% initialise
input_parameters;
Model_real = update_model(true_model_type);
iter = (usr_par.whichFrq-1) * (change_freq_every) + it +1;
% iter = usr_par.cumulative_iter + 1;

disp ' ';
disp(['---- generating output for iteration ', num2str(iter)]);
disp ' ';


% % from usr_par
output_path = usr_par.output_path;
sObsPerFreq = usr_par.sObsPerFreq;
Model_bg    = usr_par.Model_bg;
% Kg          = usr_par.Kg(iter);
% Kseis       = usr_par.Kseis(iter);
% K_total     = usr_par.K_total(iter);
 InvProps    = usr_par.InvProps;
% if isfield(usr_par, 'Model_real')
%     Model_real = usr_par.Model_real;
% end
% 
% 
%% preparation

% domain
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

% convert model to usable format
[Model] = map_m_to_parameters(m, usr_par);
K_abs   = map_gradm_to_gradparameters(gm, usr_par);

%% inversion output
% InvProps.misfit(iter) = jm;
% InvProps.misfitseis(iter) = NaN;
% InvProps.misfitgrav(iter) = NaN;

% % useful output
% if strcmp(use_grav,'no')
%     Kg{iter}=NaN;
% end
% if exist('Model_real', 'var')
%     InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model, Model_real);
% else
%     InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model);
% end


%% plot stuff

% model
fig_mod = plot_model_diff(Model,Model_bg,param_plot);
titel = [project_name,': model diff of iter ', num2str(iter), ' and bg model'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = [output_path,'/iter',num2str(iter,'%03d'),'.model-diff.',param_plot,'.png'];
% print(fig_mod,'-dpng', '-r0', figname);
print(fig_mod,'-dpng','-r400',figname);
close(fig_mod);
clearvars('fig_mod');

% load obs information belonging to the current ModRandString
ModFolder = [output_path,'/fwd_temp/',ModRandString,'/'];
load([ModFolder,'iter-rec.mat']);

% gravity difference
fig_grav_comp = plot_gravity_quivers(usr_par.rec_g, g_recIter, usr_par.g_obs, ...
                X, Z, Model.rho);
figname = [output_path,'/iter',num2str(iter,'%03d'),'.gravity_difference.png'];
titel = ['gravity diff of model - real model (iter ', num2str(iter), ')'];
mtit(fig_grav_comp, titel, 'xoff', 0.001, 'yoff', 0.00001);
% print(fig_grav_comp, '-dpng', '-r0', figname);
print(fig_grav_comp, '-dpng', '-r400', figname);
close(fig_grav_comp); clearvars('fig_grav_comp');


% seismograms
% disp 'plotting seis?'
% if strcmp(use_seis, 'yesseis')
%     for isrc = 1:nsrc
%         vel = sEventRecIter(isrc).vel;
%         vobs = usr_par.sEventObs(isrc).vel;
%         t   = sEventRecIter(isrc).t;
%         recs = 1:length(vel);
%         % determine how many seismograms are actually plotted
%         if length(vel) > 8; recs = [2:2:length(vel)]; end
%         % actual plotting
%         fig_seisdif = plot_seismogram_difference(vel, vobs, t, recs);
%         figname = [output_path,'/iter.',num2str(iter,'%03d'),'.seisdif.src',num2str(isrc),'.png'];
%         print(fig_seisdif,'-dpng','-r400',figname); close(fig_seisdif);
%     end
% end

% 
% % gravity kernel (if applicable)
% if strcmp(use_grav,'yes')
%     fig_Kg  = plot_gravity_kernel(Kg);
%     figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernel_grav.rho.png'];
%     titel = [project_name, '- Gravity kernel for iter ',num2str(iter)];
%     mtit(fig_Kg,titel, 'xoff', 0.001, 'yoff', 0.00001);
%     print(fig_Kg,'-dpng','-r400',figname);
%     close(fig_Kg);
%     clearvars fig_Kg;
% end
% 
% % seismic kernels
% % absolute rho-mu-lambda
% fig_knl = plot_kernels(Kseis, 'rhomulambda',Model, 'total', 'own', 99.95);
% titel = [project_name,' - iter ',num2str(iter), ' seismic kernels (abs rho-mu-lambda)'];
% mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
% figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels.abs.rho-mu-lambda.png'];
% print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
% % absolute rho-vs-vp
% fig_knl = plot_kernels(Kseis, 'rhovsvp',Model, 'total', 'own', 99.95);
% titel = [project_name, ' - iter ',num2str(iter), ' seismic kernels (abs rho-vs-vp)'];
% mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
% figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels.abs.rho-vs-vp.png'];
% print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
% 

% total kernel in (relative) rhomulambda
K_rel = calculate_relative_kernels(K_abs, Model_bg);
fig_knl = plot_kernels(K_rel, 'rhomulambda', Model, 'total', 'same', 99.95);
titel = [project_name,' - iter ',num2str(iter), ' TOTAL kernels (rel rho-mu-lambda)'];
mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels-total.rel.rho-mu-lambda.png'];
% print(fig_knl,'-dpng','-r0', figname); close(fig_knl);
print(fig_knl,'-dpng','-r400',figname); close(fig_knl);

% fig_knl = plot_kernels(K_total, 'rhovsvp',Model, 'total', 'own', 99.95);
% titel = [project_name,' - iter ',num2str(iter), ' TOTAL kernels (abs rho-vs-vp)'];
% mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
% figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels-total.abs.rho-vs-vp.png'];
% print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
% 
% % density kernel buildup:  seis + grav = total
% if strcmp(use_grav, 'yes')
%     
%     switch parametrisation
%         case 'rhomulambda'
%             Krho.seis = filter_kernels(Kseis.rho.total,smoothgwid);
%             Krho.grav = filter_kernels(Kg,smoothgwid);
%             Krho.together = filter_kernels(K_total.rho.total,smoothgwid);
%         case 'rhovsvp'
%             Kseis2 = change_parametrisation_kernels('rhomulambda',parametrisation,Kseis,Model_bg);
%             Ktot2  = change_parametrisation_kernels('rhomulambda',parametrisation,K_total,Model_bg);
%             Krho.seis = filter_kernels(Kseis2.rho2.total,smoothgwid);
%             Krho.grav = filter_kernels(Kg,smoothgwid);
%             Krho.together = filter_kernels(Ktot2.rho2.total,smoothgwid);
%         otherwise
%             error('the parametrisation in which kernels are added was unknown');
%     end
%     
%     fig_Krho = plot_model(Krho);
%     maks = prctile(abs([Krho.seis(:); Krho.grav(:);Krho.together(:)]),99.5);
%     for ii = 2:2:6
%         fig_Krho.Children(ii).CLim = [-maks maks];
%     end
%     titel = [project_name,' - buildup of density kernel - iter ',num2str(iter)];
%     mtit(fig_Krho,titel, 'xoff', 0.001, 'yoff', 0.04);
%     figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernel-rho-buildup.png'];
%     print(fig_Krho,'-dpng','-r400',figname); close(fig_Krho);
%     
% end
% 
% 
% 
% 
% if (iter > 1)
%         % inversion results with inversion landscape plot
%         fig_inv2 = plot_inversion_development_landscapeshape(InvProps, iter);
%         figname = [output_path,'/inversion_development.',project_name,'.misfit-landscape'];
%         print(fig_inv2,'-dpng','-r400',[figname,'.png']);
%         print(fig_inv2,'-depsc','-r400',[figname,'.eps']);
%         close(fig_inv2)
%         
%         fig_invres = plot_inversion_result(InvProps, iter);
%         titel = [project_name,': inversion results'];
%         mtit(fig_invres,titel, 'xoff', 0.0000001, 'yoff', 0.03);
%         figname = [output_path,'/inversion_result.',project_name];
%         print(fig_invres,'-dpng','-r400',[figname,'.png']);
%         print(fig_invres,'-depsc','-r400',[figname,'.eps']);
%         close(fig_invres)
% end


%% save stuff to usr_par

% iter specific info
usr_par.Model(iter) = Model;
usr_par.K_abs(iter) = K_abs;
usr_par.misfit(iter) = jm;
usr_par.cumulative_iter = iter;

%% plot inversion result 

% set previous to NaN if it doesn't exist
if iter == 1
    usr_par.previous = NaN;
end

% set real model to NaN if it doens't exist
if ~exist('Model_real', 'var')
    Model_real = NaN;
end

% load 'current iter' variables from file
ModFolder = [output_path,'/fwd_temp/',ModRandString,'/'];
load([ModFolder,'currentIter.misfits.mat']);
load([ModFolder,'currentIter.kernels.mat']);
blips = fieldnames(currentKnls);
currentIter = currentMisfits;
for ii = 1:numel(blips)
    currentIter.(blips{ii}) = currentKnls.(blips{ii});
end

% calculate some inversion output numbers: kernel magnitudes, model
% diffferences, angles between kernels, and then plot that.
InvProps = calc_inversion_output_InvTlbx(iter, currentIter, usr_par.previous, InvProps, usr_par.K_abs, usr_par.Model, jm, gm, Model_real);
if iter > 1
    fig_invres = plot_inversion_result_InvTlbx(InvProps, iter);
    titel = [project_name,' - iter ',num2str(iter), ' inversion result'];
    mtit(fig_invres,titel, 'xoff', 0.001, 'yoff', 0.04);
    figname = [output_path,'/inversion_result.png'];
%     print(fig_invres,'-dpng','-r0', figname); close(fig_invres);
    print(fig_invres,'-dpng','-r400',figname); close(fig_invres);
end

usr_par.InvProps = InvProps;


%% save output to file

% if ~exist([output_path,'/iter',num2str(iter,'%03d'),'.all-vars.mat'],'file')
disp 'saving all current variables...'
clearvars('figname', 'savename', 'fig_seisdif', 'fig_mod', ...
    'filenm_old', 'filenm_new', 'fig_knl');
savename = [output_path,'/iter',num2str(iter,'%03d'),'.all-vars.mat'];

excludedVars = {'sObsPerFreq', 'currentIter', ...
                'currentKnls', 'K_rel', 'K_abs', 'gm', ...
                'Model_real', 'Model_bg', 'Model', 'm', ...
                'X', 'Z'};

save_except(savename, excludedVars{1:end});
% end

    
%% move all output files to the actual output path
blips = dir([output_path,'*iter.current*']);
for ii = 1:numel(blips)
    bestand = blips(ii).name;
    oldfile = [output_path,bestand];
    newfile = [output_path,strrep(blips(ii).name,'iter.current',['iter',num2str(iter,'%03d')])];
    movefile(oldfile,newfile);
end




%% prepare usr_par for new information from the next iter.
% !!! NOTE IT DOES NOT WORK LIKE THIS IF A NEW FREQ IS INITIATED!!
%     because then, previous is overwritten, with information of a model
%     that is otherwise never used. Will have to think about that, but it
%     will cause WRONG results in calculating e.g. the angles between prev
%     & current model just after the frequency has changed...
%             -- Nienke Blom, 7 August 2015
usr_par.previous = currentIter;
% rmfield(usr_par, 'current');

%% throw away temp folder

TempDir = [output_path,'/fwd_temp'];
if (exist(TempDir, 'dir'))
    rmdir(TempDir, 's');
end


%% new iteration


disp ' ';
disp '=====================================';
disp(['==== starting iteration ',num2str(iter+1)]);
disp '=====================================';
disp ' ';

end
