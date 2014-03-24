% travel time kernel calculation routine

% run forward
cd ../code/
[v_rec_x,v_rec_y,v_rec_z,t,rec_x,rec_z,test]=run_forward;


cd ../tools/
% check what the seismograms look like
v_rec = [v_rec_x; v_rec_y; v_rec_z;];
plot_seismograms(v_rec,t,'velocity');

% make adjoint sources
adstf = make_all_adjoint_sources(v_rec_x,v_rec_y,v_rec_z,t);

% check the adjoint source time functions
plot_vrec_to_adjointstf(t,v_rec_x,squeeze(adstf(1,:,:)));
plot_vrec_to_adjointstf(t,v_rec_z,adstf(2,:,:));
plot_vrec_to_adjointstf(t,v_rec_y,adstf(3,:,:));

% run adjoint 
cd ../code/
K=run_adjoint('traveltime');

Lx=;     % model extension in x-direction [m]
Lz=;     % model extension in z-direction [m]
nx=;     % grid points in x-direction
nz=;     % grid points in z-direction
stf_PSV=;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
[mu,rho,lambda]=define_material_parameters(nx,nz,11);
set_figure_properties;

cd ../tools/
plot_kernels_rho_mu_lambda;

calculate_other_kernels