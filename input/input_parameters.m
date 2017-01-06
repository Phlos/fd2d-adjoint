%==========================================================================
% project name (all file names will be changed accordingly)
%==========================================================================

project_name='Systematic.test-061.anoms-10-pct';

%==========================================================================
% inversion properties
%==========================================================================

% adjoint_source_path='../input/sources/adjoint/';

% apply hard constraints?
apply_hc = 'no';   % 'yes' or 'no'
% hard constraints
axrot = 'x';     % 'x' or 'z' at the moment.

% use gravity?
use_grav = 'no'; % 'yes' or 'no'
which_grav = 'g_vector'; % 'g_vector' or 'g_potential'
use_seis = 'yesseis'; % 'yesseis' or 'noseis'

% what misfit functional are we using
misfit_type = 'waveform_difference'; % 'waveform_difference' or 'cc_time_shift
% inversion parametrisation
parametrisation = 'rhomulambda';   % 'rhovsvp' or 'rhomulambda', maybe later 'rhomukappa'
param_plot = 'rhovsvp';

% fix velocities?
fix_velocities = 'no'; % 'yes' or 'no'

% normalise misfits:
normalise_misfits = 'byfirstmisfit'; % 'byfirstmisfit' or 'div_by_obs' or 'no'
                                  % 'byfirstmisfit':
                                  % normalises both the s and g misfits by 
                                  % their first value, so that they're same
                                  % order of magnitude
                                  % 'div_by_obs':
                                  % divides L2 norm by obs values so that
                                  % you get percent wise misfit


% initial step length;
% stepInit = 3.5e14;    % good for circular configuration
% stepInit = 5e15;        % good for circular src and rec @ top of domain
% stepInit = 5e14;        % good for src @ bottom, rec @ top of domain (2-12-2014)
% stepInit = 1e-1;        % kernels normalised by 1st misfit size. (20-11-2014)
% stepInit = 1e-8;        % normalised misfit, TINY rho anomaly (model 101) (24-11-2014)
% ---------------- changed stf adstf: now divided by surface of cell to
%                  make it a spatial delta function
% stepInit = 5e9;         % KMP solved dd ~25 feb 2015
% stepInit = 7e7;         % PREM + 1% vs anomalies
% stepInit = 1e4;         % PREM + 1% rho2 anomalies
% stepInit = 1e8;         % tt and wavef inv: truemod = 100, starting = 1; (21-3-2015)
% stepInit = 5e6;         % low freq (0.01 Hz) PREM + 1000 kg/m3 (23-3-2015)
stepInit = 0.02;        % L-BFGS with kernels corrected (July 2015)
% stepInit = 0.0025;      % L-BFGS, seis only, rho-vs-vp

%- smoothing properties
% % smoothing (= filtering) seismograms before adstf
% max_freq = 0.2; % Hz

% smoothing kernels?
smoothing  = 'yessmooth'; % 'yessmooth' or 'nosmooth'
smoothgwid = 2; % width of the gaussian in the smoothing filter (pixels)
                % used to be 9 w/ conv2 
                
% zero out the bottom 5 rows of the kernel:
zero_bottom_rows = 'yeszerobottom'; % 'yeszerobottom' or 'nozerobottom'
                
% store forward wavefield every .. timesteps
store_fw_every = 10; 

%==========================================================================
% set basic simulation parameters
%==========================================================================

wave_propagation_type='PSV';   % can be 'PSV' or 'SH' or 'both'

Lx=6000e3;     % model extension in x-direction [m]
Lz=2890e3;     % model extension in z-direction [m] ! PREM: don't exceed 2891

% nx=901;     % grid points in x-direction --> PREM: dx = 6.67 km
% nz=430;     % grid points in z-direction --> PREM: dz = 6.72 km
% nx = 241;   % PREM: dx = 25 km
% nz = 116;   % PREM: dz = 25 km
% nx = 601;   % PREM: dx = 10 km
% % nz = 290;   % PREM: dz = 10 km
% nx = 301;   % PREM: dx = 20 km
% nz = 145;   % PREM: dz = 20 km
nx = 430;   % PREM: dx = ~14 km  (13.99)
nz = 207;   % PREM: dz = ~14 km  (14.03)
% nx = 1201;
% nz = 581;

% The necesssary time step (in order to obtain a stable model run) may vary
% according to the chosen gridding. 
% dt=1.0;      % time step [s] - PREM model, dx=dz=25km
% dt = 0.8;    % time step [s] - dx=dz= 20 km
dt = 0.6;     % time step [s] - dx=dz = 14 km
% dt=0.1;       % time step [s] 0.5 explodes, 0.4 suffices @PREM dx=dz=10km
tmax = 1200;    % length of run [s] -- 1200 should be enough for ScS (=935 s) (PcP = 510)
nt = ceil(tmax/dt); % number of iterations
nt=store_fw_every*round(nt/store_fw_every);

order=4;    % finite-difference order (2 or 4) (2 is not recommended)

%==========================================================================
% model type
%==========================================================================

% starting model from matfile
use_matfile_startingmodel = 'no';
% starting_model = '../output/Model_0.03Hz.mat';
starting_model = '';

bg_model_type = 50;     % PREM
true_model_type = 87;   % PREM + LM and UM solid circles (all params)
model_type=50; %(start) % PREM + LM/UM solid circles in vp, vs only

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
% 50= PREM background values model (will be selected at true height above cmb)
% 60 = PREM bg + random rho2 AND vs AND vp, 1% anomalies
% 61 = PREM bg + ONLY random vs & vp anomalies, 1%
% 82 = small solid non-overlapping hard-edged circles in rho, vs, and vp
% 83 = same as 82, but only rho circles
% 84 = same as 82, but only vs, vp circles
% 85 = PREM + LM and UM solid circles (all params)
% 86 = PREM + LM and UM solid circles -- only vs and vp
% 100= layered: left = high rho0, right = low rho0
% 101= homogeneous model (Tromp like) with tiny rho anomaly
% 102= Ring shaped model (Evangelos): vp=5000, vs=3000, rho=2600 | outside: 5000,1,2600

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)



%==========================================================================
% sources -- positions
%==========================================================================

nsrc = 8;
src_depth = 50e3; % source depth beneath surface (m)

%- line of sources that have two sources in the same place - only possible w/ even nsrc
pos_x = (1: 1: nsrc/2) * (Lx/(nsrc/2+1));
for ii = 1:nsrc
    if mod(ii,2)==0 % even
        src_info(ii).loc_x = pos_x(ii/2);
    else % odd
        src_info(ii).loc_x = pos_x((ii-1)/2 + 1);
    end
    src_x(ii) = src_info(ii).loc_x;
        
    src_info(ii).loc_z = (Lz - src_depth); % sources at 50 km depth
    src_z(ii) = src_info(ii).loc_z;
end


% %- line of sources near top/bot of the domain
% src_x= (1: 1: nsrc) * (Lx/(nsrc+1));
% % dz = 1/16 * Lz;
% dz = 50e3;
% src_z=ones(size(src_x)) * (Lz - 2*dz); % -2*dz necessary as a result of b.c.)
% % src_z=ones(size(src_x)) * (2*dz); % at the bottom, 2*dz from the bottom ivm absbound
% for ii = 1:nsrc
%     src_info(ii).loc_x = src_x(ii);
%     src_info(ii).loc_z = src_z(ii);
% end
clearvars src_depth pos_x;

%==========================================================================
% sources -- source-time functions
%==========================================================================

% with this loop, all sources are the same w/ the same polarisation (7-5-2015)
for ii = 1:nsrc
    src_info(ii).stf_type = 'delta_bp';   % 'ricker' or 'delta_bp' (twice butterworth bandpassed 
                                         % delta function)
    % needed for 'ricker'
    src_info(ii).tauw_0  = 2.628;         % seconds
    src_info(ii).tauw    = 10*4.0;        % source duration, seconds
    src_info(ii).tee_0   = 10*2.5;        % source start time, seconds
    
    % needed for 'delta_bp'    
    src_info(ii).f_min=0.006667;          % minimum stf frequency [Hz]
    src_info(ii).f_max=1.0;               % maximum stf frequency [Hz]

    if mod(ii,2)==0; % even
        src_info(ii).stf_PSV = [1 0]; % [x z] --> S waves radiate up/down
    else % odd
        src_info(ii).stf_PSV = [0 1]; % [x z] --> P waves radiate up/down
    end
                        % direction of the source-time-function in P-SV wave 
                        % propagation. The final stf will be normalised
                        % such that its original amplitude is preserved.
end
source_amplitude = 1e9;                 
                    
%- source filtering - 8 frequency bands increasing by a factor 1.25 each time
filter_stf_with_freqlist = true;
butterworth_npoles = 5;
f_minlist = [0.00667 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667];
f_maxlist = [0.00667 0.00833 0.01042 0.01302 0.01628 0.02035 0.02543 0.03179];

% how many iterations with the same source?
change_freq_every = 20;          % how many iterations with the same freq?



%==========================================================================
% receiver positions
%==========================================================================

%- a line of receivers just below the top boundary
nrec = 16;
% nrec = 1;
rec_x= (1: 1: nrec) * (Lx/(nrec+1));
dz = Lz/(nz-1);
rec_z=ones(size(rec_x)) * (Lz-2*dz); % -2*dz necessary as a result of b.c.)


%==========================================================================
% gravity measurement positions
%==========================================================================

%- a line of gravity receivers above the domain
rec_height = 20e3; % [m]
nrec_g = 50;
rec_g.x= (0: 1: nrec_g) * (Lx/nrec_g);
rec_g.z=ones(size(rec_g.x)) * (Lz + rec_height);

%==========================================================================
% simulation mode
%==========================================================================

simulation_mode='forward';

% 'forward'                 regular forward simulation
% 'forward_green'           forward simulation where Fourier transform of the Greens function is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present
% 'noise_source_kernel'     compute sensitivity kernel for the noise source power-spectral density distribution

%==========================================================================
% absorbing boundaries
%==========================================================================

width = 500.0e3;        % width of the boundary layer in m 

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=0;   % absorb waves on the top boundary
absorb_bottom=0;% absorb waves on the bottom boundary

%==========================================================================
% plotting
%==========================================================================

% plot every 'plot every'th image (otherwise computationally rather heavy)
plot_every=nt*2; % value larger than nt, so that no plotting takes place
% plot_every = 40;

plot_forward_frames='PSV';   % 'X-Y-Z' or 'X-Y' or 'PSV-SH' or 'PSV' 
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
make_movie='no';                                   % 'yes' or 'no'
make_movie_adj='no';                               % 'yes' or 'no'
movie_file=['../output/',project_name,'.forward'];        % output file name
movie_file_adj=['../output/',project_name,'.adjoint'];
