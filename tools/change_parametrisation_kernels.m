function Kout = change_parametrisation_kernels(inparam, outparam, Kin, model1, varargin)

% function that changes the parametrisation of sensitivity kernels in an
% elastic domain. Given an input parametrisation with three fields (e.g. 
% rho mu lambda), it calculates the corresponding values of the fields in 
% the output parametrisation.
% NOTE: only 'total' subkernels are converted!
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

%% preparation

% reparametrising absolute or relative kernels?
absrel = checkargs(varargin(:));

% model in rho-vs-vp to be able to perform the reparametrisations
model = change_parametrisation('rhomulambda','rhovsvp',model1);
    
  
%% actual reparametrisations:

switch outparam
    
    case 'rhovsvp'
        
        if isfield(Kin, 'rho2') && isfield(Kin, 'vs2') && isfield(Kin, 'vp2')
            warning('your infields already has the parametrisation you want');
            Kout = Kin;
        end
        
        
        switch inparam
            case 'rhomulambda'
                
                if strcmp(absrel,'rel')
                    Kin.rho.total = Kin.rho.total ./ model1.rho;
                    Kin.mu.total  = Kin.mu.total ./ model1.mu;
                    Kin.lambda.total  = Kin.lambda.total ./ model1.lambda;
                end
                
%                 fig_knlin = plot_kernels(Kin);
%                 titel = ['input kernel in rhomulambda'];
%                 mtit(fig_knlin, titel);
%                 pause(1)
                
                Kout.rho2.total = Kin.rho.total ...
                    + model.vs.^2 .* Kin.mu.total ...
                    + (model.vp .^2 - 2.*model.vs .^2) ...
                    .* Kin.lambda.total;
                Kout.vs2.total  = 2*model.rho .* model.vs .* Kin.mu.total ...
                    - 4*model.rho .* model.vs .* Kin.lambda.total;
                Kout.vp2.total  = 2*model.rho .* model.vp .* Kin.lambda.total;
%                 fig_knlout = plot_kernels(Kout);
%                 titel = ['output kernel in rhovsvp'];
%                 mtit(fig_knlout, titel);
%                 pause(15)
            otherwise
                warning('Unexpected output parametrisation.');
        end

        if strcmp(absrel,'rel')
            Kout.rho2.total = Kout.rho2.total .* model.rho;
            Kout.vs2.total  = Kout.vs2.total .* model.vs;
            Kout.vp2.total  = Kout.vp2.total .* model.vp;
        end
        

    % outparam = rhomulambda
    case 'rhomulambda'
        
        if isfield(Kin, 'mu') && isfield(Kin, 'lambda')
%             warning('your infields already has the parametrisation you want');
            Kout = Kin;
        else
            switch inparam
                case 'rhovsvp'

                    if strcmp(absrel,'rel')
                        Kin.rho2.total = Kin.rho2.total .* model.rho;
                        Kin.vs2.total  = Kin.vs2.total .* model.vs;
                        Kin.vp2.total  = Kin.vp2.total .* model.vp;
                    end
                
%                 fig_knlin = plot_kernels(Kin);
%                 titel = ['input kernel in rhovsvp'];
%                 mtit(fig_knlin, titel);
%                 pause(15);
                    Kout.rho.total = Kin.rho2.total ...
                                     - 0.5 * model.vs ./ model.rho .* Kin.vs2.total ...
                                     - 0.5 * model.vp ./ model.rho .* Kin.vp2.total;
                    Kout.mu.total = 0.5 ./ (model.rho .* model.vs) .* Kin.vs2.total ...
                                    + 1 ./ (model.rho .* model.vp) .* Kin.vp2.total;
                    Kout.lambda.total = 0.5 ./ (model.rho .* model.vp) .* Kin.vp2.total;
%                 fig_knlout = plot_kernels(Kout);
%                 titel = ['output kernel in rhomulambda'];
%                 mtit(fig_knlout, titel);
%                 pause(15);
                
                otherwise
                    warning('Unexpected input parametrisation.');
            end
        
            
        end
        
        if strcmp(absrel,'rel')
            Kout.rho.total = Kout.rho.total ./ model1.rho;
            Kout.mu.total  = Kout.mu.total ./ model1.mu;
            Kout.lambda.total  = Kout.lambda.total ./ model1.lambda;
        end

    otherwise
        warning('Unexpected output parametrisation.');
end

end

function absrel = checkargs(args)


narg = length(args);

if narg == 0
%     disp 'no absrel, we assume absolute kernels'
    absrel = 'abs';
elseif narg == 1
    absrel = args{1};
    disp(['chn_par_knl: given the provided input, we''re using ',absrel,' kernels'])
else
    disp 'unsupported number of input argumetns to change_parametrisation_kernels';
end

end