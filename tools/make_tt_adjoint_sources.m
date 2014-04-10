function [adstf] = make_tt_adjoint_sources(v_rec,t)

% wrapper to make the adjoint sources for travel time sensitivity kernels

adjoint_source_component = input('which component of P-SV do you read? [ x / z ] ', 's');

if strcmp(adjoint_source_component,'x')
    disp 'makin'' X -- Z is zero'
    [adstf_x, misfit_x]=make_adjoint_sources(v_rec.x,zeros(size(v_rec.x)),t,'velocity','cc_time_shift','_1');
    adstf_z = zeros(size(adstf_x));
elseif strcmp(adjoint_source_component,'z')
    disp 'makin'' Z -- X is zero'
    [adstf_z, misfit_z]=make_adjoint_sources(v_rec.z,zeros(size(v_rec.z)),t,'velocity','cc_time_shift','_3');
    adstf_x = zeros(size(adstf_z));
end

adstf_SH = input('Do you want to calculate an SH adjoint? [yes / no] ', 's');
if (strcmp(adstf_SH,'y') || strcmp(adstf_SH,'yes'))
    [adstf_y, misfit_y]=make_adjoint_sources(v_rec.y,zeros(size(v_rec.y)),t,'velocity','cc_time_shift','_2');
elseif (strcmp(adstf_SH,'n') || strcmp(adstf_SH,'no'))
    error('I did not understand your input')
end

adstf(1,:,:) = adstf_x;
adstf(2,:,:) = adstf_y;
adstf(3,:,:) = adstf_z;