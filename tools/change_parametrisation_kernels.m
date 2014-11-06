function Kout = change_parametrisation_kernels(inparam, outparam, Kin, model1)

% function that changes the parametrisation of sensitivity kernels in an
% elastic domain. Given an input parametrisation with three fields (e.g. 
% rho mu lambda), it calculates the corresponding values of the fields in 
% the output parametrisation.
% 
% NOTE: For now, only rhovsvp <--> rhomulambda works, but I can easily expand.
%
% INPUT:
% inparam:  parametrisation of the input kernels: 'rhovsvp' or 'rhomulambda'
% outparam: parametrisation of the output kernels, idem
% in:       input kernels, structured as (at least):
%           in.param1.total, in.param2.total, in.param3.total
% model:    the current model, structured as:
%           model.rho, .mu, .lambda, .vs, .vp
%
% OUTPUT:
% out:      output kernels, structured as:
%           out.param1.total, out.param2.total, out.param3.total
%
% -- Nienke Blom, 27-10-2014

model = change_parametrisation('rhomulambda','rhovsvp',model1);

switch outparam
    case 'rhovsvp'
        if isfield(Kin, 'rho2') && isfield(Kin, 'vs2') && isfield(Kin, 'vp2')
            warning('your infields already has the parametrisation you want');
%             Kout = Kin;
        end
            switch inparam
                case 'rhomulambda'
                    Kout.rho2.total = Kin.rho.total ...  
                                      + model.vs.^2 .* Kin.mu.total ...
                                      + (model.vp .^2 - 2.*model.vs .^2) ...
                                        .* Kin.lambda.total;
                    Kout.vs2.total  = 2*model.rho .* model.vs .* Kin.mu.total ...
                                      - 4*model.rho .* model.vs .* Kin.lambda.total;
                    Kout.vp2.total  = 2*model.rho .* model.vp .* Kin.lambda.total;
                otherwise
                    warning('Unexpected output parametrisation.');
            end
%         end
    case 'rhomulambda'
        
        if isfield(Kin, 'mu') && isfield(Kin, 'lambda')
            warning('your infields already has the parametrisation you want');
%             Kout = Kin;
            
        end
        
            switch inparam
                case 'rhovsvp'
                    Kout.rho.total = Kin.rho2.total ...
                                     - 0.5 * model.vs ./ model.rho .* Kin.vs2.total ...
                                     - 0.5 * model.vp ./ model.rho .* Kin.vp2.total;
                    Kout.mu.total = 0.5 ./ (model.rho .* model.vs) .* Kin.vs2.total ...
                                    + 1 ./ (model.rho .* model.vp) .* Kin.vp2.total;
                    Kout.lambda.total = 0.5 ./ (model.rho .* model.vp) .* Kin.vp2.total;
                otherwise
                    warning('Unexpected output parametrisation.');
            end
%         end
    otherwise
        warning('Unexpected input parametrisation.');
end

end