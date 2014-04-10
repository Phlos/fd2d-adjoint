function [adstf, misfit] = make_all_adjoint_sources(v_rec,v_obs,t,measurement)

% wrapper to make the adjoint sources for waveform inversion
if strcmp(measurement,'waveform_difference')
    size(v_rec.y)
    [adstf_x, misfit.x] = make_adjoint_sources(v_rec.x,v_obs.x, ...
                          t,'displacement',measurement,'_1');
    [adstf_z, misfit.z] = make_adjoint_sources(v_rec.z,v_obs.z, ...
                          t,'displacement',measurement,'_3');
    
    adstf_SH = input('Do you want to calculate an SH adjoint? [yes / no] ', 's');
    if (strcmp(adstf_SH,'y') || strcmp(adstf_SH,'yes'))
        [adstf_y, misfit.y]=make_adjoint_sources(v_rec.y,v_obs.y, ...
                            t,'displacement',measurement,'_2');
    elseif (strcmp(adstf_SH,'n') || strcmp(adstf_SH,'no'))
        adstf_y = zeros(size(v_rec.y));
        %     misfit_y= 0;
    else
        error('I did not understand your input')
    end
    
    adstf(1,:,:) = adstf_x;
    adstf(2,:,:) = adstf_y;
    adstf(3,:,:) = adstf_z;
    
else
    error('I did not understand the adjoint measurement you want to make.')
end