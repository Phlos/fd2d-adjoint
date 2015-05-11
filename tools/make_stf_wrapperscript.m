function [stf, t] = make_stf_wrapperscript(f_min_fromlist, f_max_fromlist)
% makes the iter dependent = frequency dependent source time function

%% prepare
input_parameters;

sfe = store_fw_every;
nt=sfe*round(nt/sfe);
t=0:dt:dt*(nt-1);

%% calculate
% f_max_fromiterlist
switch stf_type
        case {'delta_bp', 'heaviside_bp'}
            stf = make_source_time_function(t,stf_type,f_min,f_max);
        case 'ricker'
            error('we have not yet obtained a frequency dependent ricker')
            stf = make_source_time_function(t,stf_type,tauw_0, tauw, tee_0);
            % bandpass the ricker wavelet!
            
end

end