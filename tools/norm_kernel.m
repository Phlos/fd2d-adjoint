function K_normd = norm_kernel(kernel, misfit_1, normalise_how)

% Normalises a kernel by misfit_1 -- if applicable
% - normalise_how = 'byfirstmisfit': divide misfit by 1st misfit
% - normalise_how = 'no':            misfit_normd = misfit

if isstruct(kernel)
    params = fieldnames(kernel);
    for m = 1:length(params)
        if isstruct(kernel.(params{m}))
            comps = fieldnames(kernel.(params{m}));
            for l = 1:length(comps);
                K_in = kernel.(params{m}).(comps{l}); % shorthand
                K_out = actual_normalisation(K_in, misfit_1, normalise_how);
                K_normd.(params{m}).(comps{l}) = K_out; % shorthand
            end
        else
            error('the input kernel to norm_kernel is neither seis nor grav?!')
        end
    end
    
elseif isnumeric(kernel)
    K_normd = actual_normalisation(kernel, misfit_1, normalise_how);
else
    error('the input kernel to norm_kernel is neither seis nor grav?!')
end



end

function K_normd = actual_normalisation(kernel, misfit_1, normalise_how)

switch normalise_how
    case 'byfirstmisfit'
        if misfit_1 == 0
            warning('your initial misfit is zero, cannot divide by that!')
            K_normd = kernel;
        else
            K_normd = kernel ./ misfit_1;
        end
    case 'no'
        K_normd = kernel;
end

end