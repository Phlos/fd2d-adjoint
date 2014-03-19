% travel time kernel calculation routine

cd ../code/
[v_rec_x,v_rec_y,v_rec_z,t,rec_x,rec_z]=run_forward;
cd ../tools/
misfit_x=make_adjoint_sources(v_rec_x,zeros(size(v_rec_x)),t,'velocity','cc_time_shift','_1');
misfit_y=make_adjoint_sources(v_rec_y,zeros(size(v_rec_y)),t,'velocity','cc_time_shift','_2');
misfit_z=make_adjoint_sources(v_rec_z,zeros(size(v_rec_z)),t,'velocity','cc_time_shift','_3');
cd ../code/
K=run_adjoint('traveltime');
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
cd ../tools/
plot_kernels_rho_mu_lambda;