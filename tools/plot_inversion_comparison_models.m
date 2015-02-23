% script to print start & end models in the right way

% getting minmax values
minmax.rho(1) = 2.53814e+03;
minmax.rho(2) = 2.66220e+03;
minmax.vs(1) = 3.12332e+03;
minmax.vs(2) =  3.27070e+03;
minmax.vp(1) = 5.68583e+03;
minmax.vp(2) = 5.90973e+03;

% plotting the models
fig.invseis.end = plot_model_minmax(invseis_Model(80),minmax,'rhovsvp');
fig.invseis.start = plot_model_minmax(invseis_Model(1),minmax,'rhovsvp');
fig.invseisgrav.start = plot_model_minmax(invseisgrav_Model(1),minmax,'rhovsvp');
fig.invseisgrav.end = plot_model_minmax(invseisgrav_Model(100),minmax,'rhovsvp');
fig.invseisgravhcx.start = plot_model_minmax(invseisgravhcx_Model(1),minmax,'rhovsvp');
fig.invseisgravhcx.end = plot_model_minmax(invseisgravhcx_Model(100),minmax,'rhovsvp');
fig.inv_truemodel = plot_model_minmax(Model_real,'rhovsvp');

% printing the models
set(fig.invseis.start,'PaperOrientation','landscape')
set(fig.invseis.start,'PaperUnits','normalized')
set(fig.invseis.start,'PaperPosition', [0 0 1 1])
set(fig.invseis.start,'Renderer','painters');
figname = '../output/inv-seis.Model-start.pdf';
print(fig.invseis.start,'-dpdf','-r200',figname);


set(fig.invseis.end,'PaperOrientation','landscape')
set(fig.invseis.end,'PaperUnits','normalized')
set(fig.invseis.end,'PaperPosition', [0 0 1 1])
set(fig.invseis.end,'Renderer','painters');
figname = '../output/inv-seis.Model-end.pdf';
print(fig.invseis.end,'-dpdf','-r200',figname);

set(fig.invseisgrav.start,'PaperOrientation','landscape')
set(fig.invseisgrav.start,'PaperUnits','normalized')
set(fig.invseisgrav.start,'PaperPosition', [0 0 1 1])
set(fig.invseisgrav.start,'Renderer','painters');
figname = '../output/inv-seis-grav.Model-start.pdf';
print(fig.invseisgrav.start,'-dpdf','-r400',figname);

set(fig.invseisgrav.end,'PaperOrientation','landscape')
set(fig.invseisgrav.end,'PaperUnits','normalized')
set(fig.invseisgrav.end,'PaperPosition', [0 0 1 1])
set(fig.invseisgrav.end,'Renderer','painters');
figname = '../output/inv-seis-grav.Model-end.pdf';
print(fig.invseisgrav.end,'-dpdf','-r400',figname);

set(fig.invseisgravhcx.start,'PaperOrientation','landscape')
set(fig.invseisgravhcx.start,'PaperUnits','normalized')
set(fig.invseisgravhcx.start,'PaperPosition', [0 0 1 1])
set(fig.invseisgravhcx.start,'Renderer','painters');
figname = '../output/inv-seis-grav-hcx.Model-start.pdf';
print(fig.invseisgravhcx.start,'-dpdf','-r400',figname);

set(fig.invseisgravhcx.end,'PaperOrientation','landscape')
set(fig.invseisgravhcx.end,'PaperUnits','normalized')
set(fig.invseisgravhcx.end,'PaperPosition', [0 0 1 1])
set(fig.invseisgravhcx.end,'Renderer','painters');
figname = '../output/inv-seis-grav-hcx.Model-end.pdf';
print(fig.invseisgravhcx.end,'-dpdf','-r400',figname);

% true model
set(fig.inv_truemodel,'PaperOrientation','landscape')
set(fig.inv_truemodel,'PaperUnits','normalized')
set(fig.inv_truemodel,'PaperPosition', [0 0 1 1])
set(fig.inv_truemodel,'Renderer','painters');
figname = '../output/inv-all.Model-true.pdf';
print(fig.inv_truemodel,'-dpdf','-r400',figname);

% set(fig.invseis.start,'Renderer','painters');
% figname = '../output/inv-seis.Model-start.eps';
% print(fig.invseis.start,'-depsc','-r400',figname);
% 
% set(fig.invseis.end,'Renderer','painters');
% figname = '../output/inv-seis.Model-end.eps';
% print(fig.invseis.end,'-depsc','-r400',figname);
% 
% set(fig.invseisgrav.start,'Renderer','painters');
% figname = '../output/inv-seis-grav.Model-start.eps';
% print(fig.invseisgrav.start,'-depsc','-r400',figname);
% 
% set(fig.invseisgrav.end,'Renderer','painters');
% figname = '../output/inv-seis-grav.Model-end.eps';
% print(fig.invseisgrav.end,'-depsc','-r400',figname);
% 
% set(fig.invseisgravhcx.start,'Renderer','painters');
% figname = '../output/inv-seis-grav-hcx.Model-start.eps';
% print(fig.invseisgravhcx.start,'-depsc','-r400',figname);
% 
% set(fig.invseisgravhcx.end,'Renderer','painters');
% figname = '../output/inv-seis-grav-hcx.Model-end.eps';
% print(fig.invseisgravhcx.end,'-depsc','-r400',figname);


% figname = '../output/inv-seis.Model-end-r200.eps';
% print(fig.invseis.end,'-depsc','-r200',figname);
% figname = '../output/inv-seis.Model-end-r200.pdf';
% print(fig.invseis.end,'-dpdf','-r200',figname);
% set(fig.invseis.end,'PaperUnits','normalized')
% figure(fig.invseis.end)
% figure(figtest)
% figure(fig.invseis.end)
% 
% figtest = plot_model_minmax(invseis_Model(80),minmax,'rhovsvp');
% set(figtest,'PaperUnits','normalized')
% set(figtest,'PaperPosition', [0 0 1 1])
% set(figtest,'Renderer','painters')
% set(figtest,'PaperOrientation','landscape')
% print(figtest,'-dpdf','-r200','../output/figtest.pdf');