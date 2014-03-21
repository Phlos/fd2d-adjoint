% travel time kernel calculation routine

cd ../code/
[v_rec_x,v_rec_y,v_rec_z,t,rec_x,rec_z]=run_forward;
cd ../tools/
[adstf_x, misfit_x]=make_adjoint_sources(v_rec_x,zeros(size(v_rec_x)),t,'velocity','cc_time_shift','_1');
[adstf_y, misfit_y]=make_adjoint_sources(v_rec_y,zeros(size(v_rec_y)),t,'velocity','cc_time_shift','_2');
% [adstf_z, misfit_z]=make_adjoint_sources(v_rec_z,zeros(size(v_rec_z)),t,'velocity','cc_time_shift','_3');

% check the adjoint source time functions
plot_vrec_to_adjointstf(t,v_rec_x,adstf_x);
plot_vrec_to_adjointstf(t,v_rec_z,adstf_z);
plot_vrec_to_adjointstf(t,v_rec_y,adstf_y);

cd ../code/
K=run_adjoint('traveltime');
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
cd ../tools/
plot_kernels_rho_mu_lambda;
[mu,rho,lambda]=define_material_parameters(nx,nz,11);
calculate_other_kernels