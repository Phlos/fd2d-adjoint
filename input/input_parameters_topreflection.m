%==========================================================================
% project name (all file names will be changed accordingly)
%==========================================================================

project_name='top-reflection';

%==========================================================================
% path where seismic sources are located
%==========================================================================

adjoint_source_path='../input/sources/adjoint/';
% adjoint_source_component='x';   % 'x' or 'z' -- the P-SV source direction
                                % the sensitivity kernel is calculated
                                % based on the component given here. 

%==========================================================================
% set basic simulation parameters
%==========================================================================

wave_propagation_type='PSV';   % can be 'PSV' or 'SH' or 'both'

Lx=2.20e5;     % model extension in x-direction [m]
Lz=1.1e5;     % model extension in z-direction [m]

nx=660;     % grid points in x-direction
nz=330;     % grid points in z-direction

% The necesssary time step (in order to obtain a stable model run) may vary
% according to the chosen gridding. 
% dt=0.33;     % time step [s] fine for SH in a grid   dx=dz=3.34e3 m
% dt=0.1;      % time step [s] fine for P-SV in a grid dx=dz=3.34e3 m
dt=0.01;      % time step [s]
nt=5000;     % number of iterations

order=4;    % finite-difference order (2 or 4)

%==========================================================================
% model type
%==========================================================================

model_type=11;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient medium
% 6=vertical gradient medium with localised density perturbation
% 11= homogeneous with values from Tromp et al 2005
% 100= layered: left = high velocity, right = low velocity (any difference with model 3???)

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)

%==========================================================================
% source-time function
%==========================================================================

stf_type = 'ricker';    % 'ricker' or 'delta_bp' (twice butterworth bandpassed 
                        % delta function)
% needed for 'ricker'
tauw_0  = 2.628;      % seconds
tauw    = 4.0;        % source duration, seconds
tee_0   = 2.5;        % source start time, seconds

% needed for 'delta_bp'    
f_min=0.2;     % minimum frequency [Hz]
f_max=1.00;     % maximum frequency [Hz]

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

src_x=[0.6e5];
src_z=[0.7e5];

%==========================================================================
% receiver positions
%==========================================================================

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
rec_x=[1.6e5];
rec_z=[0.7e5];

%- a large number of receivers in a closed rectangular configuration
%rec_x=[50.0  50.0  50.0  50.0  50.0   50.0    70.0  90.0 110.0 130.0   70.0  90.0 110.0 130.0  150.0 150.0 150.0 150.0 150.0  150.0];
%rec_z=[70.0  90.0 110.0 130.0 150.0  170.0    70.0  70.0  70.0  70.0  170.0 170.0 170.0 170.0   70.0  90.0 110.0 130.0 150.0  170.0];


%==========================================================================
% absorbing boundaries
%==========================================================================

width=25000.0;     % width of the boundary layer in m

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=0;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary

%==========================================================================
% plotting
%==========================================================================

% plot every 'plot every'th image (otherwise computationally rather heavy)
plot_every=25;

plot_forward_frames='X-Y';   % 'X-Y-Z' or 'X-Y' or 'PSV-SH' or 'PSV' 
                             % which frames should be plotted in the forward calculation

%==========================================================================
% output: movies, matfiles, etc.
%==========================================================================

%- screenshots of wave propagation

snapshotfile = ['../output/',project_name];
savetimes = [0 5 7 10 15 20 25 30 35 40 45];
% savetimes = [];

%- matfiles ---

save_u_fw = 'no';       % 'yes' or 'no' -- save the u_forward matfile
save_v_fw = 'no';       % 'yes' or 'no' -- save the v_forward matfile


%- movies -----

make_movie='yes';                                   % 'yes' or 'no'
make_movie_adj='yes';                               % 'yes' or 'no'
movie_file=['../output/',project_name,'_forward'];        % output file name
movie_file_adj=['../output/',project_name,'_adjoint'];