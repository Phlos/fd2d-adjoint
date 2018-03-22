
% function stf = make_source_time_function(t,stf_type,simulation_mode,f_min,f_max,tauw_0,tauw,tee_0)
% function stf = make_source_time_function(t, src_amp, dx, dz, stf_PSV, stf_type, varargin)
function stf = make_source_time_function(src_number, varargin)

% compute source-time function with directional components x,y,z
%
% SYNTAX:
% stf = make_source_time_function(t,'delta_bp',f_min, f_max)
% stf = make_source_time_function(t,'ricker', tauw_0, tauw, tee_0)
%
% INPUT:
% src_number = which source do we want to produce?
%
% OUTPUT:
% stf = stf.x
%       stf.y
%       stf.z
%       struct with all components of the wavefield (Y only used for SH, X
%       and Z only for P-SV wave propagation)

% adapted on 20-3-2015, Nienke Blom
% adapted on 22-3-2018, Nienke Blom (make it a 3 component stf)

input_parameters;

[~,~,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

t = 0:dt:dt*(nt-1);
dt = t(2) - t(1);

% retrieve source information struct
src_inf = src_info(src_number);
stf_type = src_inf.stf_type;
stf_PSV = src_inf.stf_PSV;
% for stf_type 'delta_bp'
f_min = src_inf.f_min;
f_max = src_inf.f_max;
% for stf_type 'ricker'
tauw_0 = src_inf.tauw_0;
tauw   = src_inf.tauw;
tee_0  = src_inf.tee_0;
% [f_min, f_max, tauw_0, tauw, tee_0] = checkargs(stf_type, varargin(:));


switch stf_type
    case 'delta'
        stfn = zeros(1, length(t));
        stfn(1) = 1 / dt;
    
    case 'delta_bp'
% if (strcmp(stf_type,'delta_bp'))
    
        stfn=zeros(1,length(t));
        stfn(1)=1;
        stfn=butterworth_lp(stfn,t,5,f_max,'silent');
        stfn=butterworth_hp(stfn,t,3,f_min,'silent');
        
        % normalise amplitude
        stfn = stfn ./ max(abs(stfn));
    
% elseif (strcmp(stf_type,'ricker'))
    case 'ricker'
        
        alfa = 2 * tauw_0 / tauw;
        stfn = ( -2*alfa^3 / pi) * (t-tee_0) .* exp( -alfa^2 * (t-tee_0).^2);
        
        % normalise amplitude
        stfn = stfn ./ max(abs(stfn));
        
    case 'heaviside_bp'
        
%         stf=zeros(1,length(t)); stf(round(length(t)/2):length(t)) = 1;
        stfn=ones(1,length(t)); 
        stfn=butterworth_lp(stfn,t,3,f_max,'silent');
        stfn=butterworth_hp(stfn,t,3,f_min,'silent');
        
    otherwise
        error('stf_type not recognised');
        
end

% adapt source time function amplitude to input
stfn = source_amplitude * stfn;

% % should make this into a plot all srces x and z (and y)
% fig_stf = plot_source_time_function(t,stfn);
% close(fig_stf);

% prefactor = so that the stf is a spatial delta function (integral 1)
% correct for point-source-ness of source: divide by dx dz to
% make independent of grid size
prefac = 1.0 / dx / dz;
stfn = prefac * stfn;

%- insert source time function into x y z with proper magnitudes
stf.x = stfn .* stf_PSV(1)./norm(stf_PSV);  % x direction
stf.y = stfn;                               % y direction
stf.z = stfn .* stf_PSV(2)./norm(stf_PSV);  % z direction

% if strcmp(simulation_mode,'forward_green')
%     
%     %- To compute the Fourier transform of the Greens function, we need
%     %- a Heaviside function here. This is because the code computes the
%     %- velocity field Greens function and not the displacement field
%     %- Greens function. So, we integrate by integrating the source time
%     %- function.
%     
%     stf=1.0e9*ones(1,length(t));
% end




end
