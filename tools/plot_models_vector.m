%- plot all figures nicely

% inversion development
fig_inv = plot_inversion_development(InvProps, ...
                           Model(niter), Model(1), Model_real, middle);
figname = ['../output/inversion_development.',project_name,'.pdf'];
set(fig_inv,'Renderer','painters')
print(fig_inv,'-dpdf','-r400',figname);
close(fig_inv)

% real model
if(exist('Model_real','var'))
    fig_mod = plot_model_diff(Model_real, Model_start, param_plot);
    set(fig_mod,'Renderer','painters')
    figname = ['../output/obs.',param_plot,'.eps'];
    print(fig_mod,'-depsc','-r400',figname);
    close(fig_mod)
else
    warning('It seems like the real model doesn''t exist...')
end

% model(1)
fig_mod = plot_model_diff(Model(1),Model_start,param_plot);
set(fig_mod,'Renderer','painters');
titel = [project_name,': model of iter 1'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter001.model.',param_plot,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
close(fig_mod);

% model(end)
fig_mod = plot_model_diff(Model(niter),Model_start,param_plot);
set(fig_mod,'Renderer','painters');
titel = [project_name,': final model'];
mtit(fig_mod, titel, 'xoff', 0.001, 'yoff', -0.05);
figname = ['../output/iter',num2str(niter),'.model.',param_plot,'.eps'];
print(fig_mod,'-depsc','-r400',figname);
close(fig_mod);