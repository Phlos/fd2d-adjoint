function sEventInfo = prepare_stf()

% prepares a source-time function w/ locations & src direction & everything
%
% SYNTAX:
% sEventInfo = prepare_stf();
%
% INPUT:
% information gotten from input_parameters
%
% OUTPUT:
% sEventInfo:   struct like sEventInfo{isrc}.loc_x
%                                           .loc_z
%                                           .stf  .x
%                                                 .y
%                                                 .z
%                                           .t
%
% -- N.A. Blom, 13 May 2015

%% prep

% get from input:
%   - list of source locations
%   - list of source directions
%   - list of source time function functions
% -> this is all in src_info(1...nsrc)
input_parameters;
[~,~,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

sfe = store_fw_every;
nt=sfe*round(nt/sfe);
t = 0:dt:dt*(nt-1);

for ii = 1:length(src_info)
    
    % copy relevant info from src_info to sources
    sEventInfo(ii).loc_x = src_info(ii).loc_x;
    sEventInfo(ii).loc_z = src_info(ii).loc_z;
    
    % prepare info for stf calculation
    stf_type = src_info(ii).stf_type;
    f_min = src_info(ii).f_min; f_max = src_info(ii).f_max;
    tauw_0 = src_info(ii).tauw_0; tauw = src_info(ii).tauw; tee_0 = src_info(ii).tee_0;
    stf_PSV = src_info(ii).stf_PSV;
    
    
    %% stf calculation
    % make the actual source time functions for each source
    switch stf_type
        case {'delta'}
            stfn = make_source_time_function(t, stf_type);
        case {'delta_bp', 'heaviside_bp'}
            stfn = make_source_time_function(t,stf_type,f_min,f_max);
        case 'ricker'
            stfn = make_source_time_function(t,stf_type,tauw_0, tauw, tee_0);
    end
    stfn = source_amplitude * stfn;
    
    % should make this into a plot all srces x and z (and y)
    fig_stf = plot_source_time_function(t,stfn);
    

    % prefactor = so that the stf is a spatial delta function (integral 1)
    % correct for point-source-ness of source: divide by dx dz to
    % make independent of grid size
    prefac = 1.0 / dx / dz;
    stfn = prefac * stfn;
    
    %- insert source time function into x y z with proper magnitudes
    sEventInfo(ii).stf.x = stfn .* stf_PSV(1)./norm(stf_PSV);  % x direction
    sEventInfo(ii).stf.y = stfn;                             % y direction
    sEventInfo(ii).stf.z = stfn .* stf_PSV(2)./norm(stf_PSV);  % z direction
    
    sEventInfo(ii).t = t;

    close(fig_stf);
end

end