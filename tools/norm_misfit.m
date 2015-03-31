function [misfit_normd, div_by_first, div_by_obs] = ...
    norm_misfit(misfit,normalise_how, misfit_1, obs)

% Normalises a misfit by misfit_1 -- if applicable
% - normalise_how = 'byfirstmisfit': divide misfit by 1st misfit
% - normalise_how = 'no':            misfit_normd = misfit

%% the division factors
% calculate sum_obs
sum_obs = 0;
if iscell(obs)
    for ii = 1:length(obs)
        comp = fieldnames(obs{ii});
        for icomp = 1:length(comp)
            sum_obs = sum_obs + sum(obs{ii}.(comp{icomp}) .^2);
        end
    end
else
    sum_obs = sum(obs.x(:) .^2 ) + sum(obs.z(:) .^2);
end

if misfit_1 == 0;
%  disp 'dividing by 1st misfit and by sum_obs1 == 0
%             warning('your initial misfit is zero, cannot divide by that!')
    div_by_first = misfit;
else
    div_by_first = misfit / misfit_1;
end
if sum_obs == 0;
    eor('sum(v_obs^2) is zero, cannot divide by that! AND shouldn''t be!')
    div_by_obs = misfit;
else
    div_by_obs = misfit / sum_obs;
end

switch normalise_how
%masaving the correct one as misfit_normdhow
    case 'byfirstmisfit'
        misfit_normd = div_by_first;
% %         misfit_1 = normdiv;
%         if misfit_1 == 0
% %             warning('your initial misfit is zero, cannot divide by that!')
%             misfit_normd = misfit;
%         else
%             misfit_normd = misfit / misfit_1;
%         end
    case 'div_by_obs'
        misfit_normd = div_by_obs;
%         sumvobs = normdiv;
%         if sumvobs == 0
%             error('sum(v_obs^2) is zero, cannot divide by that! AND shouldn''t be!')
%             misfit_normd = misfit;
%         else
%             misfit_normd = misfit / misfit_1;
%         end
    case 'no'
        misfit_normd = misfit;
end

end

