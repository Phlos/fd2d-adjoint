function outfields = change_parametrisation(inparam, outparam, infields)

% function that changes the parametrisation of an elastic domain. Given an
% input parametrisation with three fields (e.g. rho mu lambda), it
% calculates the corresponding values of the fields in the output
% parametrisation.
% For now, only rhovsvp <--> rhomulambda works, but I can easily expand.
%
% -- Nienke Blom, 18-8-2014

switch inparam
    case 'rhomulambda'
        switch outparam
            case 'rhovsvp'
                outfields.vs = abs(sqrt(infields.mu ./ infields.rho));
                outfields.vp = abs(sqrt((infields.lambda + 2*infields.mu) ...
                               ./ infields.rho));
                outfields.rho = infields.rho;
            otherwise
                warning('Unexpected output parametrisation.');
        end
    case 'rhovsvp'
        switch outparam
            case 'rhomulambda'
                outfields.mu  = abs(infields.vs .^2 .* infields.rho);
                outfields.lambda = abs(infields.rho2 .* ...
                                   ( infields.vp.^2 - 2* infields.vs.^2));
                outfields.rho = infields.rho;
            otherwise
                warning('Unexpected output parametrisation.');
        end
    otherwise
        warning('Unexpected input parametrisation.');
end
        

end