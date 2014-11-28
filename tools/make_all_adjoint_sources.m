function [adstf, misfit] = make_all_adjoint_sources(v_rec,v_obs,t,measurement,mode,varargin)

% v_rec         recorded velocity field at all receivers
% v_obs         observed velocity field at all receivers
% t             time axis
% measurement   'waveform_difference' or 'cc_time_shift' -- the first for
%               full waveform inversion, the second for traveltime kernels.
% mode          'auto' or 'manual': do you want to pick your seismograms by
%               hand or should the computer use & taper the whole thing?'

% get scalingfactor from possible varargin
scalingfactor = checkargs(varargin(:));

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
    adstf_SH = 'no';
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
    
    
    % calculate the total misfit
    comp = fieldnames(misfit);
    misfit.total = 0;
    for component = comp'
        misfit.total = misfit.total + misfit.(component{1});
    end
    % determine scalingfactor if the misfit has to be scaled by itself
    % (iter 1)
    if isnan(scalingfactor)
        scalingfactor = misfit.total;
    end
    misfit.normd = misfit.total / scalingfactor;
    
else
    error('I did not understand the adjoint measurement you want to make.')
end

end

function sf = checkargs(args)

% checking input to make_all_adjoint_sources
% no arguments:     sf = 1 (no scaling)
% args{1} = NaN:    sf = NaN, so that the misfit will be scaled by itself
% args{1} = number: sf = args{1}, misfit will be scaled by that number.

nargs = length(args);

if (nargs == 0)
%     disp 'make_all_adjoint_sources: no scaling factor supplied'
    sf = 1;
elseif (nargs == 1 && isnumeric(args{1}))
    sf = args{1};
%     disp(['make_all_adjoint_sources: scaling factor = ', num2str(sf)]);
    if isnan(sf)
        disp ' --> misfit will be scaled by itself'
    end
else
    error('too many or wrong input arguments to make_all_adjoint_sources')
end


end