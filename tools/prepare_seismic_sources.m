function seismic_sources = prepare_seismic_sources()

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

sfe = store_fw_every;
nt=sfe*round(nt/sfe);
t = 0:dt:dt*(nt-1);

for ii = 1:length(src_info)
    
    % copy relevant info from src_info to sources
    seismic_sources(ii).loc_x = src_info(ii).loc_x;
    seismic_sources(ii).loc_z = src_info(ii).loc_z;
    
    %% stf calculation
    % make the actual source time functions for each source
    stf = make_source_time_function(ii); % three component
    
    seismic_sources(ii).stf = stf;
    
    seismic_sources(ii).t = t;
end

end