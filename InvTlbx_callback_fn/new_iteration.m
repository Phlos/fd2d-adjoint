function [usr_par] = new_iteration( m, it, usr_par) 
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
iter = it+1;

% % from usr_par
output_path = usr_par.output_path;
sObsPerFreq = usr_par.sObsPerFreq;
Model_bg    = usr_par.Model_bg;
% Kg          = usr_par.Kg(iter);
% Kseis       = usr_par.Kseis(iter);
% K_total     = usr_par.K_total(iter);
% InvProps    = usr_par.InvProps;
% if isfield(usr_par, 'Model_real')
%     Model_real = usr_par.Model_real;
% end
% 
% 
%% load model in usable format
[Model] = map_m_to_parameters(m, usr_par);
% 
% %% inversion output
% 
% % useful output
% if strcmp(use_grav,'no')
%     Kg{iter}=NaN;
% end
% if exist('Model_real', 'var')
%     InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model, Model_real);
% else
%     InvProps = calc_inversion_output(iter, InvProps, K_total, Kg, Kseis, Model);
% end
% 
% 
%% plot stuff

% model
fig_mod = plot_model_diff(Model,Model_bg,param_plot);
titel = [project_name,': model diff of iter ', num2str(iter), ' and bg model'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = [output_path,'/iter',num2str(iter,'%03d'),'.model-diff.',param_plot,'.png'];
print(fig_mod,'-dpng','-r400',figname);
close(fig_mod);
clearvars('fig_mod');


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
% % total kernel in (absolute) rhovsvp
% fig_knl = plot_kernels(K_total, 'rhomulambda',Model, 'total', 'own', 99.95);
% titel = [project_name,' - iter ',num2str(iter), ' TOTAL kernels (abs rho-mu-lambda)'];
% mtit(fig_knl,titel, 'xoff', 0.001, 'yoff', 0.04);
% figname = [output_path,'/iter',num2str(iter,'%03d'),'.kernels-total.abs.rho-mu-lambda.png'];
% print(fig_knl,'-dpng','-r400',figname); close(fig_knl);
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
%     
% %% save all output files to the actual output path
% blips = dir('../output/*iter.current*');
% for ii = 1:numel(blips)
%     bestand = blips(ii).name;
%     oldfile = ['../output/',bestand];
%     newfile = [output_path,strrep(blips(ii).name,'iter.current',['iter',num2str(iter,'%03d')])];
%     movefile(oldfile,newfile);
% end

%% initialise new frequency if necessary



cfe = change_freq_every;
whichFrq = floor((iter-1)/cfe)+1;
if whichFrq > length(sObsPerFreq)
    whichFrq = length(sObsPerFreq);
end

sEventInfo = usr_par.sObsPerFreq(whichFrq).sEventInfo;
sEventObs   = usr_par.sObsPerFreq(whichFrq).sEventObs;


%% save stuff to usr_par

usr_par.Model(iter) = Model;
% usr_par.misfit(iter) = usr_par.misfit_total;
usr_par.sEventInfo  = sEventInfo;
usr_par.sEventObs   = sEventObs;


end
