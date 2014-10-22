% plot all the diff plots

% load compare model
load('output/ParamRhoVsVp_RhoAnom.hc_rhovsvp/ParamRhoVsVp_RhoAnom.hc_rhovsvp.all-vars.mat','Model_hcrhovsvp')

% rho at tipping vs. rho at non-tipping point
fig_notip_tip = plot_model_diff(Model_nontipping(7), Model_tippingpt(7),'rhovsvp');
mtit(fig_notip_tip,'Diff plot: inversion with rho anomaly at non MoI tipping point - anomaly at tipping point (both iter 7).');
figname = ['../output/compare.model.diff-notipping-tipping.rhovsvp.iter7.png'];
print(fig_notip_tip, '-dpng', '-r400', figname);


% real model vs. non-tipping point inversion result (iter 7)
fig_real_notip = plot_model_diff(Model_real, Model_nontipping(7), 'rhovsvp');
mtit(fig_real_notip,'Diff plot: real model - iter 7 of inversion with rho anomaly at non MoI tipping point.');
figname = ['../output/compare.model.diff-real-notipping.rhovsvp.iter7.png'];
print(fig_real_notip, '-dpng', '-r400', figname);

