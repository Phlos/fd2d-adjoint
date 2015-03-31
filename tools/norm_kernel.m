function K_normd = norm_kernel(kernel, normalise_how, misfit_1, obs)

% Normalises a kernel by misfit_1 -- if applicable
% - normalise_how = 'byfirstmisfit': divide misfit by 1st misfit
% - normalise_how = 'no':            misfit_normd = misfit

sum_obs = 0;
if iscell(obs)
    for ii = 1:length(obs)
        comp = fieldnames(obs{ii});
        for icomp = 1:length(comp)
            sum_obs = sum_obs + sum(obs{ii}.(comp{icomp}) .^2);
        end
    end
else
    sum_obs = sum(obs.mag(:) .^2);
%     sum_obs = sum(obs.x(:) .^2) + sum(obs.z(:) .^2);
end

switch normalise_how
    case 'byfirstmisfit'
        div_by = misfit_1;
    case 'div_by_obs'
        div_by = sum_obs;
    case 'no'
        div_by = 1;
    otherwise
        error('Panic! No logical normalise_how!')
end
        

if isstruct(kernel)
    params = fieldnames(kernel);
    for m = 1:length(params)
        if isstruct(kernel.(params{m}))
            comps = fieldnames(kernel.(params{m}));
            for l = 1:length(comps);
                K_in = kernel.(params{m}).(comps{l}); % shorthand
                K_out = actual_normalisation(K_in, div_by);
                K_normd.(params{m}).(comps{l}) = K_out; % shorthand
            end
        else
            error('the input kernel to norm_kernel is neither seis nor grav?!')
        end
    end
    
elseif isnumeric(kernel)
    K_normd = actual_normalisation(kernel, div_by);
else
    error('the input kernel to norm_kernel is neither seis nor grav?!')
end



end

function K_normd = actual_normalisation(kernel, div_by)

        if div_by == 0
            warning('your divide_by is zero, cannot divide by that!')
            K_normd = kernel;
        else
            K_normd = kernel ./ div_by;
        end

end