%==========================================================================
% project name (all file names will be changed accordingly)
%==========================================================================

project_name='Evangelos';

%==========================================================================
% inversion properties
%==========================================================================

adjoint_source_path='../input/sources/adjoint/';

misfit_type = 'cc_time_shift'; % 'waveform_difference' or 'cc_time_shift'

% apply hard constraints?
apply_hc = 'yes';   % 'yes' or 'no'
% hard constraints
axrot = 'x';     % 'x' or 'z' at the moment.

% use gravity?
use_grav = 'yes'; % 'yes' or 'no'

% inversion parametrisation
parametrisation = 'rhovsvp';   % 'rhovsvp' or 'rhomulambda', maybe later 'rhomukappa'

% normalise misfits:
normalise_misfits = 'byfirstmisfit'; % 'byfirstmisfit' or 'no' % normalises both the s and g misfits by their first value, so that they're same magnitudes

% initial step length;
% stepInit = 3.5e14;    % good for circular configuration
% stepInit = 5e15;        % good for circular src and rec @ top of domain
% stepInit = 5e14;        % good for src @ bottom, rec @ top of domain (2-12-2014)
stepInit = 1e-1;        % kernels normalised by 1st misfit size. (20-11-2014)
% stepInit = 1e-8;        % normalised misfit, TINY rho anomaly (model 101) (24-11-2014)

% smoothing properties
% smoothnp = 15;  % size of the smoothing filter
smoothgwid = 5; % width of the gaussian in the smoothing filter (pixels)
                % used to be 9 w/ conv2 




%==========================================================================
% set basic simulation parameters
%==========================================================================

wave_propagation_type='SH';   % can be 'PSV' or 'SH' or 'both'

Lx=1.50e5;     % model extension in x-direction [m]
Lz=1.50e5;     % model extension in z-direction [m]

nx=450;     % grid points in x-direction
nz=450;     % grid points in z-direction

% The necesssary time step (in order to obtain a stable model run) may vary
% according to the chosen gridding. 
% dt=0.33;     % time step [s] fine for SH in a grid   dx=dz=3.34e3 m
% dt=0.1;      % time step [s] fine for P-SV in a grid dx=dz=3.34e3 m
dt=0.01;      % time step [s]
nt=10000;      % number of iterations

order=4;    % finite-difference order (2 or 4) (2 is not recommended)

%==========================================================================
% model type
%==========================================================================

model_type=10;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient mu medium
% 6=vertical gradient mu medium with localised density perturbation
% 10= homogeneous with values from Tromp et al 2005
% 11= like 10, but with a positive rectangular mu anomaly in the centre
% 12= like 11, but now it's a density anomaly.
% 14= gaussian rho anomaly in the centre.
% 15= gaussian mu anomaly in the centre.
% 17= gaussian rho anomaly off-centre.
% 18= gaussian mu anomaly off-centre.
% 21= gaussian off-central rho2 anomaly (rho2 = rho in rho-vs-vp)
% 31= five 'rand' positive rho2 anomalies (rho2 = rho in rho-vs-vp)
% 41= ten 'rand' rho2 anomalies: 5 pos 5 neg (rho2 = rho in rho-vs-vp)
% 100= layered: left = high velocity, right = low velocity (any difference with model 3???)

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)


% % set the background values for plot_model for all models based on Tromp
% trompmodels = [10:39];
% if any(model_type==trompmodels)
%     disp 'real model based on Tromp'
%     middle.rml = [2600 2.66e10    2.42e10];
%     middle.rvv = [2600 3198.55736 5797.87759];
% end


%==========================================================================
% source-time function
%==========================================================================

stf_type = 'ricker';    % 'ricker' or 'delta_bp' (twice butterworth bandpassed 
                        % delta function)
% needed for 'ricker'
tauw_0  = 2.628;      % seconds
tauw    = 4.0;        % source duration, seconds
tee_0   = 2.5;        % source start time, seconds
% tauw = 4.0;
% tee_0 = 4.0;

% needed for 'delta_bp'    
f_min=0.2;          % minimum frequency [Hz]
f_max=1.00;         % maximum frequency [Hz]

stf_PSV = [1 0];    % [x z]
                    % direction of the source-time-function in P-SV wave 
                    % propagation. The final stf will be normalised with
                    % such that its original amplitude is preserved.

%==========================================================================
% simulation mode
%==========================================================================

simulation_mode='forward';

% 'forward'                 regular forward simulation
% 'forward_green'           forward simulation where Fourier transform of the Greens function is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present
% 'noise_source_kernel'     compute sensitivity kernel for the noise source power-spectral density distribution

%==========================================================================
% source positions
%==========================================================================

% %- line of sources at the bottom of the domain -- use with absbound bottom?
% nsrc = 8;
% % nsrc = 1;
% src_x= (1: 1: nsrc) * (Lx/(nsrc+1));
% dz = 1/16 * Lz;
% src_z=ones(size(src_x)) * (0 + 2*dz); % -2*dz necessary as a result of b.c.)

% %- 'random' source positions, 8x
% src_x = Lx *   [ 0.2769    0.0462    0.0971    0.8235    0.6948    0.3171    0.9502    0.0344];
% src_z = Lz * (1-[0.4387    0.3816    0.7655    0.7952    0.1869    0.8981    0.4456    0.6463]);

% %- circle around the middle of the domain
% centre=[Lx/2 Lz/2];
% numrec=8;
% circlesize = 0.25*min(Lx,Lz);
% dphi = 2*pi/(numrec);
% 
% src_x=zeros(1,numrec);
% src_z=zeros(1,numrec);
% 
% n=1;
% for phi = -pi : dphi : pi-dphi ;
%     src_x(n)=centre(1) + circlesize*cos(phi);
%     src_z(n)=centre(2) + circlesize*sin(phi);
%     n=n+1;
% end

% % a single source
src_x = Lx/8;
src_z = Lz/2;
% src_x=[0.6e5];
% src_z=[0.7e5];

%==========================================================================
% receiver positions
%==========================================================================

% %- a line of receivers just below the top boundary
% nrec = 16;
% % nrec = 1;
% rec_x= (1: 1: nrec) * (Lx/(nrec+1));
% dz = Lz/(nz-1);
% rec_z=ones(size(rec_x)) * (Lz-2*dz); % -2*dz necessary as a result of b.c.)

%- a circle of receivers with the centre at centre
centre=[Lx/2 Lz/2];
numrec=8;
dphi = 2*pi / numrec;
circlesize = 0.375*min(Lx,Lz);

rec_x=zeros(1,numrec);
rec_z=zeros(1,numrec);

n=1;
for phi=-pi+dphi/2 : dphi : pi-dphi/2 ;
% for phi=-pi : dphi : pi-dphi ;
    rec_x(n)=centre(1) + circlesize*cos(phi);
    rec_z(n)=centre(2) + circlesize*sin(phi);
    n=n+1;
end

% a set of receivers in 1/3 circle around the source (hardcoded distances!)
% rec_x=zeros(1,6);
% rec_z=zeros(1,6);
% n=1;
% for phi=-pi/4:pi/10:pi/2
%     rec_x(n)=src_x(1)+2.5e5*cos(phi);
%     rec_z(n)=src_z(1)+2.5e5*sin(phi);
%     n=n+1;
% end

%- a large number of receivers in an open rectangular configuration
%rec_x=[50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0  50.0  50.0  50.0  60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 60.0  70.0  80.0  90.0  100.0 110.0 120.0 130.0];
%rec_z=[70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0 170.0 180.0 70.0 70.0 70.0 70.0 70.0  70.0  70.0  70.0  180.0 180.0 180.0 180.0 180.0 180.0 180.0 180.0];

%- just one receiver
% rec_x=[1.6e5];
% rec_z=[0.7e5];
% rec_x = 1.5e5;
% rec_z = 0.6e5;

%- a large number of receivers in a closed rectangular configuration
%rec_x=[50.0  50.0  50.0  50.0  50.0   50.0    70.0  90.0 110.0 130.0   70.0  90.0 110.0 130.0  150.0 150.0 150.0 150.0 150.0  150.0];
%rec_z=[70.0  90.0 110.0 130.0 150.0  170.0    70.0  70.0  70.0  70.0  170.0 170.0 170.0 170.0   70.0  90.0 110.0 130.0 150.0  170.0];


%==========================================================================
% gravity measurement positions
%==========================================================================

%- a line of gravity receivers above the domain
rec_height = 200; % m
nrec_g = 50;
rec_g.x= (0: 1: nrec_g) * (Lx/nrec_g);
rec_g.z=ones(size(rec_g.x)) * (Lz + rec_height);



%==========================================================================
% absorbing boundaries
%==========================================================================

width=25000.0;     % width of the boundary layer in m

absorb_left   =0;  % absorb waves on the left boundary
absorb_right  =0;  % absorb waves on the right boundary
absorb_top    =0;  % absorb waves on the top boundary
absorb_bottom =0;  % absorb waves on the bottom boundary

%==========================================================================
% plotting
%==========================================================================

% plot every 'plot every'th image (otherwise computationally rather heavy)
% plot_every=100000; % value larger than nt, so that no plotting takes place
plot_every = 10;

plot_forward_frames='SH';   % 'X-Y-Z' or 'X-Y' or 'PSV-SH' or 'PSV' 
                             % which frames should be plotted in the forward calculation
% some test about plotting the frames differently
% plot_frame.PSV='no';
% plot_frame.SH='yes';
% plot_frame.X='yes';
% plot_frame.Z='no';

%==========================================================================
% output: movies, matfiles, etc.
%==========================================================================

%- screenshots of wave propagation

snapshotfile = ['../output/',project_name];
% savetimes = [0 5 7 10 15 20 25 27];
savetimes = [];

%- matfiles ---

save_u_fw = 'no';       % 'yes' or 'no' -- save the u_forward matfile
save_v_fw = 'no';       % 'yes' or 'no' -- save the v_forward matfile


%- movies -----
make_movie='yes';                                   % 'yes' or 'no'
make_movie_adj='no';                               % 'yes' or 'no'
movie_file=['../output/',project_name,'.forward'];        % output file name
movie_file_adj=['../output/',project_name,'.adjoint'];