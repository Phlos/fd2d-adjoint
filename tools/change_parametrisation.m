function outfields = change_parametrisation(inparam, outparam, infields)

% function that changes the parametrisation of an elastic domain. Given an
% input parametrisation with three fields (e.g. rho mu lambda), it
% calculates the corresponding values of the fields in the output
% parametrisation.
% For now, only rhovsvp <--> rhomulambda works, but I can easily expand.
%
% INPUT:
% infields:     infields.param1, infields.param2, infields.param3
%
% OUTPUT:
% outfields:    outfields.param1, outfields.param2, outfields.param3
%
% -- Nienke Blom, 18-8-2014

switch outparam
    case 'rhovsvp'
        if isfield(infields, 'vs') && isfield(infields, 'vp')
            warning('your infields already has the parametrisation you want');
            outfields = infields;
        else
            switch inparam
                case 'rhomulambda'
                    outfields.rho = infields.rho;
                    outfields.vs = abs(sqrt(infields.mu ./ infields.rho));
                    outfields.vp = abs(sqrt((infields.lambda + 2*infields.mu) ...
                        ./ infields.rho));
                otherwise
                    warning('Unexpected output parametrisation.');
            end
        end
    case 'rhomulambda'
        
        if isfield(infields, 'mu') && isfield(infields, 'lambda')
            warning('your infields already has the parametrisation you want');
            outfields = infields;
            
        else
        
            switch inparam
                case 'rhovsvp'
                    outfields.rho = infields.rho;
                    outfields.mu  = abs(infields.vs .^2 .* infields.rho);
                    outfields.lambda = abs(infields.rho .* ...
                        ( infields.vp.^2 - 2* infields.vs.^2));
                otherwise
                    warning('Unexpected output parametrisation.');
            end
        end
    otherwise
        warning('Unexpected input parametrisation.');
end
        

end