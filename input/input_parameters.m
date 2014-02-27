%==========================================================================
% path where seismic sources are located
%==========================================================================

adjoint_source_path='../input/sources/adjoint/';

%==========================================================================
% set basic simulation parameters
%==========================================================================

wave_propagation_type='both';   % can be 'PSV' or 'SH' or 'both'

Lx=1.5e6;     % model extension in x-direction [m]
Lz=1.0e6;     % model extension in z-direction [m]

nx=450;     % grid points in x-direction
nz=300;     % grid points in z-direction

% The necesssary time step (in order to obtain a stable model run) may vary
% according to the chosen gridding. 
% dt=0.33;     % time step [s] fine for SH in a grid   dx=dz=3.34e3 m
% dt=0.1;      % time step [s] fine for P-SV in a grid dx=dz=3.34e3 m
dt=0.1;      % time step [s]
nt=1000;     % number of iterations (2000 reaches the receivers at dt=0.1)

order=4;    % finite-difference order (2 or 4)

%==========================================================================
% model type
%==========================================================================

model_type=1;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient medium
% 6=vertical gradient medium with localised density perturbation
% 100= layered: left = high velocity, right = low velocity (any difference with model 3???)

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)

%==========================================================================
% source-time function
%==========================================================================

f_min=0.02;     % minimum frequency [Hz]
f_max=0.10;     % maximum frequency [Hz]

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

src_x=[Lx*1/3];
src_z=[Lz*1/2];

%==========================================================================
% receiver positions
%==========================================================================

% a set of receivers in 1/3 circle around the source (hardcoded distances!)
rec_x=zeros(1,6);
rec_z=zeros(1,6);
n=1;
for phi=-pi/4:pi/10:pi/2
    rec_x(n)=src_x(1)+2.5e5*cos(phi);
    rec_z(n)=src_z(1)+2.5e5*sin(phi);
    n=n+1;
end

%- a large number of receivers in an open rectangular configuration
%rec_x=[50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0  50.0  50.0  50.0  60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 60.0  70.0  80.0  90.0  100.0 110.0 120.0 130.0];
%rec_z=[70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0 170.0 180.0 70.0 70.0 70.0 70.0 70.0  70.0  70.0  70.0  180.0 180.0 180.0 180.0 180.0 180.0 180.0 180.0];

%- just one receiver
%rec_x=[700000.0];
%rec_z=[500000.0];

%- a large number of receivers in a closed rectangular configuration
%rec_x=[50.0  50.0  50.0  50.0  50.0   50.0    70.0  90.0 110.0 130.0   70.0  90.0 110.0 130.0  150.0 150.0 150.0 150.0 150.0  150.0];
%rec_z=[70.0  90.0 110.0 130.0 150.0  170.0    70.0  70.0  70.0  70.0  170.0 170.0 170.0 170.0   70.0  90.0 110.0 130.0 150.0  170.0];


%==========================================================================
% absorbing boundaries
%==========================================================================

width=150000.0;     % width of the boundary layer in m

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=1;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary

%==========================================================================
% plotting
%==========================================================================

% plot every 'plot every'th image (otherwise computationally rather heavy)
plot_every=10;

%==========================================================================
% make wavepropagation movie
%==========================================================================

make_movie='no';                           % 'yes' or 'no'
movie_file='../output/model1_P-SV-SH';     % output file name, should be .mp4