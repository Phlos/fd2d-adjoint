%==========================================================================
% project name (all file names will be changed accordingly)
%==========================================================================

project_name='Wave-prop.S-past-vs-block';

%==========================================================================
% set basic simulation parameters
%   - model domain extent
%   - grid size in space AND time (these are related through stability
%     conditions)
%==========================================================================

wave_propagation_type='PSV';   % can be 'PSV' or 'SH' or 'both'

Lx=3000e3;     % model extension in x-direction [m]
Lz=1500e3;     % model extension in z-direction [m] ! PREM: don't exceed 2891

nx = 601;   % PREM: dx = 10 km
nz = 301;   % PREM: dz = 10 km
% these result in dx and dz (i.e. dx = Lx / (nx-1) )


% The necesssary time step (in order to obtain a stable model run) may vary
% according to the chosen gridding. 
dt = 0.1; 
%     -->  dt = 0.1 is good for visualising in dx/dz = ~10km, when the 
%          velocities are around vs=3200 m/s, vp=5800 m/s)
%     -->  @ Earth's mantle, 0.5 explodes, 0.4 suffices @PREM w/ dx=dz=10km
tmax = 500;    % length of run [s] 
%     -->  @ Earth's mantle, tmax = 1200 should be enough for ScS (=935 s)
%          (PcP = 510)
nt = ceil(tmax/dt); % number of time steps


order=4;    % finite-difference order (2 or 4) (2 is not recommended)

%==========================================================================
% model type (used in inversion, or as default if none given)
%==========================================================================

% setting the model type is normally only used in an inversion setup, but
% all the model types are given here.

% NOTE: rho2 = density parametrised with vs, vp independent
%       rho0 = density parametrised with mu, lambda independent

% 1=homogeneous with lambda = mu
% 2=homogeneous with with lambda = mu, localised rho0 perturbation
% 3=layered medium
% 4=layered with localised rho0 perturbation
% 5=vertical gradient mu medium
% 6=vertical gradient mu medium with localised density perturbation

% ==>  all models below have background values like model 10
% 10 = homogeneous with values from Tromp et al 2005 (vp=5800, vs=3200)
% -----> rho-mu-lambda parametrisation:
% 11 = block mu anomaly @ centre (+1e10)
% 12 = block rho0 anomaly @ centre (+1e3, i.e. ~40 pct) (mu, lambda const.)
% 14 = gaussian rho0 anomaly in the centre
% 15 = gaussian mu anomaly in the centre
% 17 = gaussian rho0 anomaly off-centre, low position (+1e3) (mu, lambda constant)
% 18 = gaussian mu anomaly off-centre
% -----> rho-vs-vp parametrisation:
% 13 = block rho2 anomaly @ centre (+1e3, i.e. ~40 pct) (vs, vp const.)
% 19 = block vs anomaly @ centre (+2%)
% 20 = block rho2 anomaly @ centre (+2%) (vs, vp constant)
% 22 = block rho2 anomaly @ centre (10%) (vs, vp constant)
% 23 = block vs anomaly @ centre (10%)
% 16 = gaussian central rho2 anomaly (+1e3) (vs, vp constant)
% 21 = gaussian rho2 anomaly off-centre, high (+1e3) (vs, vp constant)
% 31 = five 'random' gaussian positive rho2 anomalies (+1e3) (vs, vp constant)
% 41 = ten 'random' gaussian rho2 anomalies: 5 pos 5 neg (+/- 1e3)  (vs, vp constant)

% ==>  all models below are based on PREM -- as such, the model domain size
%      will be determined as height above CMB so really the only sensible
%      heigth is 2891 
% 50 = PREM background values model, slightly adapted
% -----> rho-mu-lambda parametrisation:
% 55 = ten 'random' gaussian rho0 anomalies: 5 pos 5 neg (+1e3) (mu, lambda const.)
% -----> rho-vs-vp parametrisation:
% 51 = ten 'random' gaussian rho2 anomalies: 5 pos 5 neg (+1e3)
% 52 = ten 'random' gaussian vs anomalies: 5 pos 5 neg (+1e3)
% 53 = ten 'random' gaussian vs anomalies: 5 pos 5 neg (+1% of max vs)
% 54 = ten 'random' gaussian rho2 anomalies: 5 pos 5 neg (+1% of max rho)
% 56 = block rho6 anomaly @ centre (+2000) (mu, lambda constant)
% 60 = ten 'random' anomalies of all three params 
       % (all +1% of their resp. max value)
       % (all parameters in different locations)
% 61 = ten 'random' anomalies of only vs & vp (1% of their resp. max values)
% 62 = ten 'random' anomalies of only vs & vp (0.5% of their resp. max values)
% 63 = ten 'random' anomalies of only vs & vp (0.9% of their resp. max values)
% 64 = ten 'random' anomalies of only vs & vp (0.95% of their resp. max values)

% ==>  hard-edged circles models (also based on PREM)
% 70 = six hard-edged circles in each parameter in a -/+/- and +/-/+ grid
% 80 = two hard-edged circles in each parameter - pos and neg
% 81 = two hard-edged circles in rho - pos and neg
% 82 = two small hard-edged circles in rho, vs, and vp
% 83 = same as 82, but only rho circles
% 84 = same as 82, but only vs, vp circles

% ==>  used for Blom et al, GJI 2017, doi:10.1093/gji/ggx076 
%      (these models are also PREM based)
% 85 = separate LM and UM solid circles, all params (+/- 1%)
% 86 = separate LM and UM solid circles, only vs and vp (+/- 1%)
% 87 = separate LM and UM solid circles, all params (+/- 10%)
% 88 = separate LM and UM solid circles, only vs and rho:
       % vs three columns + 1%, 
       % rho scaled to vs as drho/dvs = 0.4, 0.2, -0.2 .
       % This model was used for the scaling ratio tests in Blom et al.
% 89 = separate LM and UM solid circles, all params (LM: 1%, UM 10%)
% 90 = one column of LM and UM solid circles, overlapping all three params:
       % no impedance contrast, i.e. rho*vs and rho*vp stay constant over 
       % the whole model domain
       
% --> misc models
% 100= layered: left = high rho0, right = low rho0
% 101= homogeneous model (Tromp like) with tiny rho anomaly
% 102= Ring shaped model (Evangelos): vp=5000, vs=3000, rho=2600 | outside: 5000,1,2600
% 103 = PREM + raised 670 km anomaly (raised by 30 km)



%==========================================================================
% absorbing boundaries
%==========================================================================

width = 500.0e3;        % width of the boundary layer in m 

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=1;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary


%==========================================================================
% sources 
%==========================================================================

%--------------------------------------------------------------------------
% source positions 

% 1 src @ centre depth
src_info.loc_x = 1000e3;
src_info.loc_z = Lz / 2;
src_x = src_info.loc_x;
src_z = src_info.loc_z;

% nsrc = 8;
% src_depth = 50e3; % source depth beneath surface (m)
% 
% %- line of sources that have two sources in the same place 
%    -- this is only possible w/ an even nsrc (number of sources)
% pos_x = (1: 1: nsrc/2) * (Lx/(nsrc/2+1));
% for ii = 1:nsrc
%     if mod(ii,2)==0 % even
%         src_info(ii).loc_x = pos_x(ii/2);
%     else % odd
%         src_info(ii).loc_x = pos_x((ii-1)/2 + 1);
%     end
%     src_x(ii) = src_info(ii).loc_x;
%         
%     src_info(ii).loc_z = (Lz - src_depth); % sources at 50 km depth
%     src_z(ii) = src_info(ii).loc_z;
% end


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

% clearvars src_depth pos_x;

%--------------------------------------------------------------------------
% source time functions 

% with this loop, all sources are the same w/ the same polarisation (7-5-2015)
for ii = 1:length(src_info)
    src_info(ii).stf_type = 'ricker';   % 'ricker' or 'delta_bp' (twice butterworth bandpassed 
                                         % delta function)
    % needed for 'ricker'
    src_info(ii).tauw_0  = 2.628;         % seconds
    src_info(ii).tauw    = 15*4.0;        % source duration, seconds
    src_info(ii).tee_0   = 15*2.5;        % source start time, seconds
    
    % needed for 'delta_bp'    
    src_info(ii).f_min=0.006667;          % minimum stf frequency [Hz]
    src_info(ii).f_max=1.0;               % maximum stf frequency [Hz]

%     if mod(ii,2)==0; % even
%         src_info(ii).stf_PSV = [1 0]; % [x z] --> S waves radiate up/down
%     else % odd
        src_info(ii).stf_PSV = [0 1]; % [x z] --> P waves radiate up/down
%     end
                        % direction of the source-time-function in P-SV wave 
                        % propagation. The final stf will be normalised
                        % such that its original amplitude is preserved.
end
source_amplitude = 1e9;                 
                    

%==========================================================================
% receiver positions
%==========================================================================

% single receiver
% 1 receiver @ centre depth
rec_x = 2000e3;
rec_z = Lz / 2;


% %- a line of receivers just below the top boundary
% nrec = 16;
% % nrec = 1;
% rec_x= (1: 1: nrec) * (Lx/(nrec+1));
% dz = Lz/(nz-1);
% rec_z=ones(size(rec_x)) * (Lz-2*dz); % -2*dz necessary as a result of b.c.)


%==========================================================================
% gravity measurement positions
%==========================================================================

%- a line of gravity receivers above the domain
rec_height = 20e3; % [m]
nrec_g = 50;
rec_g.x= (0: 1: nrec_g) * (Lx/nrec_g);
rec_g.z=ones(size(rec_g.x)) * (Lz + rec_height);


%==========================================================================
% Plotting
%==========================================================================

%--------------------------------------------------------------------------
% wave propagation plotting

% plot every 'plot_every'th image (otherwise computationally rather heavy)
% plot_every=nt*2; % value larger than nt, so that no plotting takes place
plot_every = 25;

plot_forward_frames='PSV';   % 'X-Y-Z' or 'X-Y' or 'PSV-SH' or 'PSV' 
                             % which frames should be plotted in the forward calculation

% some test about plotting the frames differently
% plot_frame.PSV='no';
% plot_frame.SH='yes';
% plot_frame.X='yes';
% plot_frame.Z='no';


%- snapshots of wave propagation

snapshotfile = ['../output/',project_name];
% savetimes = [0 5 7 10 15 20 25 27];
savetimes = [];

%--------------------------------------------------------------------------
%- model plotting 

param_plot = 'rhomulambda'; % 'rhomulambda' or 'rhovsvp'
plot_UM_separate = false;   % only useful for PREM like models

% determine model category for plotting purposes
%  -- if the model is any of the PREM models, Z coordinate should be depth
%  beneath the surface which is set at 2891 km height... otherwise, it
%  should just be Z pointing up positive.
model_category = 'PREM_type'; % 'PREM_type' or 'other'

%- movies 
make_movie='no';                                   % 'yes' or 'no'
make_movie_adj='no';                               % 'yes' or 'no'
movie_file=['./output/',project_name,'/Velocity.forward'];        % output file name
movie_file_adj=['./output/',project_name,'/Velocity.adjoint'];
movie_label = 'S-wave past a block v_s anomaly (\rho & v_p constant)';
plot_contour_rho_anom = false;
plot_contour_vs_anom = true;


%==========================================================================
% Forward fields
%    usually stored for adjoint purposes, but can also be used for
%    snapshot and video making purposes.
%==========================================================================

% store forward wavefield every .. timesteps
store_fw_every = 10; 

% adjust the number of time steps such that storing fwd field works.
nt=store_fw_every*round(nt/store_fw_every); % rounding in case we are 
                                            % storing the fwd wavefield

%- matfiles --- (if you actually want to save to file)

save_u_fw = 'no';       % 'yes' or 'no' -- save the u_forward matfile
save_v_fw = 'no';       % 'yes' or 'no' -- save the v_forward matfile





%==========================================================================
% Inversion settings
%==========================================================================

% what misfit functional are we using
misfit_type = 'waveform_difference'; % 'waveform_difference' or 'cc_time_shift
% inversion parametrisation
parametrisation = 'rhomulambda';   % 'rhovsvp' or 'rhomulambda', maybe later 'rhomukappa'

% fix parameters?
fix_velocities = 'no'; % 'yes' or 'no'
fix_density    = 'no'; % 'yes' or 'no'

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
stepInit = 0.002;        % L-BFGS with kernels corrected (July 2015)


% smoothing kernels?
smoothing  = 'yessmooth'; % 'yessmooth' or 'nosmooth'
smoothgwid = 2; % width of the gaussian in the smoothing filter (pixels)
                % used to be 9 w/ conv2 
                
% zero out the bottom 5 rows of the kernel:
zero_bottom_rows = 'yeszerobottom'; % 'yeszerobottom' or 'nozerobottom'

% adjoint_source_path='../input/sources/adjoint/';

% apply hard constraints?
apply_hc = 'no';   % 'yes' or 'no'
% hard constraints
axrot = 'x';     % 'x' or 'z' at the moment.

% use gravity?
use_grav = 'no'; % 'yes' or 'no'
which_grav = 'g_vector'; % 'g_vector' or 'g_potential'
use_seis = 'yesseis'; % 'yesseis' or 'noseis'

%--------------------------------------------------------------------------
% MODEL TYPE: 

% starting model from matfile
use_matfile_startingmodel = 'no'; % 'yes' or 'no'
    % --> if 'yes', starting_model below should point to a model file in
    %     the correct format
% starting_model = '../output/Model_0.03Hz.mat';
starting_model = './models/Model86_perc_UM100_LM90.mat';

bg_model_type = 10;  
true_model_type = 13; 
model_type=10; %(start) 


%--------------------------------------------------------------------------
% MULTI-SCALE APPROACH -- inverting first low, then high freq bands

%- source filtering - 8 frequency bands increasing by a factor 1.25 each time
f_minlist = [0.00667]; % 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667];
f_maxlist = 0.03179;
% f_maxlist = [0.00667 0.00833 0.01042 0.01302 0.01628 0.02035 0.02543 0.03179];

% how many iterations with the same source?
change_freq_every = 1;          % how many iterations with the same freq?


%- smoothing properties
% % smoothing (= filtering) seismograms before adstf
% max_freq = 0.2; % Hz




%==========================================================================
% Deprecated? simulation mode
%==========================================================================

simulation_mode='forward';

% 'forward'                 regular forward simulation
% 'forward_green'           forward simulation where Fourier transform of the Greens function is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present
% 'noise_source_kernel'     compute sensitivity kernel for the noise source power-spectral density distribution
