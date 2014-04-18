function [adstf, misfit] = make_all_adjoint_sources(v_rec,v_obs,t,measurement,mode)

% path(path,'../input/input_parameters');
path(path,'../output/wrong_K.mu.SH_calculation/Set_006_inversion-test/')
% test3_input_parameters_testmodel;
% t=0:dt:dt*(nt-1);

% v_rec         recorded velocity field at all receivers
% v_obs         observed velocity field at all receivers
% t             time axis
% measurement   'waveform_difference' or 'cc_time_shift' -- the first for
%               full waveform inversion, the second for traveltime kernels.

% wrapper to make the adjoint sources for waveform inversion
if strcmp(measurement,'waveform_difference')
%     size(v_rec.y)
[adstf_x, misfit.x] = make_adjoint_sources(v_rec.x,v_obs.x, ...
                        t,'displacement',measurement,'_1',mode);
close(gcf);
[adstf_z, misfit.z] = make_adjoint_sources(v_rec.z,v_obs.z, ...
                        t,'displacement',measurement,'_3',mode);
close(gcf);

% adstf_SH = input('Do you want to calculate an SH adjoint? [yes / no] ', 's');
adstf_SH = 'yes';
    if (strcmp(adstf_SH,'y') || strcmp(adstf_SH,'yes'))
        [adstf_y, misfit.y]=make_adjoint_sources(v_rec.y,v_obs.y, ...
                            t,'displacement',measurement,'_2',mode);
        close(gcf);
    elseif (strcmp(adstf_SH,'n') || strcmp(adstf_SH,'no'))
        adstf_y = zeros(size(v_rec.x));
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