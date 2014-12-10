%- plot all figures nicely

% inversion development
fig_inv = plot_inversion_development(misfit, misfit_seis, misfit_g, step, modeldifn, ...
                           Model(niter), Model(1), Model_real, middle);
figname = ['../output/inversion_development.',project_name,'.pdf'];
set(fig_inv,'Renderer','painters')
print(fig_inv,'-dpdf','-r400',figname);
close(fig_inv)

% real model
if(exist('Model_real','var'))
    fig_mod = plot_model(Model_real, middle, parametrisation);
    set(fig_inv,'Renderer','painters')
    figname = ['../output/obs.',parametrisation,'.eps'];
    print(fig_mod,'-depsc','-r400',figname);
    close(fig_mod)
else
    warning('It seems like the real model doesn''t exist...')
end

% model(1)
fig_mod = plot_model(Model(1),middle,parametrisation);
set(fig_mod,'Renderer','painters');
titel = [project_name,': model of iter 1'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter001.model.',parametrisation,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
close(fig_mod);

% model(end)
fig_mod = plot_model(Model(niter),middle,parametrisation);
set(fig_mod,'Renderer','painters');
titel = [project_name,': model of iter 1'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter',num2str(niter),'.model.',parametrisation,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
close(fig_mod);