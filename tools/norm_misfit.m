function misfit_normd = norm_misfit(misfit,misfit_1,normalise_how)

% Normalises a misfit by misfit_1 -- if applicable
% - normalise_how = 'byfirstmisfit': divide misfit by 1st misfit
% - normalise_how = 'no':            misfit_normd = misfit

switch normalise_how
    case 'byfirstmisfit'
        if misfit_1 == 0
%             warning('your initial misfit is zero, cannot divide by that!')
            misfit_normd = misfit;
        else
            misfit_normd = misfit / misfit_1;
        end
    case 'no'
        misfit_normd = misfit;
end

end

