function fig_knl = plot_kernels(K, varargin)
    
% plots kernels in all formats
% 
% SYNTAX:
% fig_knl = plot_kernels(K)
% fig_knl = plot_kernels(K, parametrisation, Model)
% fig_knl = plot_kernels(K, parametrisation, Model, which_wavefields)
% fig_knl = plot_kernels(K, parametrisation, Model, which_wavefields, ...
%                        sameorownpercentile, prctiel)
% fig_knl = plot_kernels(K, parametrisation, Model, ...
%                        sameorownpercentile, prctiel)
% fig_knl = plot_kernels(K, which_wavefields);
% fig_knl = plot_kernels(K, which_wavefields, sameorownpercentile, prctiel);
% fig_knl = plot_kernels(K, sameorownpercentile, prctiel);
%
% INPUT:
% K:                    kernel (struct: K.(parameter).(wavefield)
% parametrisation:      'rhomulambda' or 'rhovsvp')
% Model:                model (struct: Model.rho, Model.mu, Model.lambda)
% which_wavefields:     'PSV', 'SH', 'total', 'both';
% sameorownpercentile:  'same' or 'own'
%                       should the kernels be plotted on the same colour
%                       scale or should it be calculated for each kernel separately?
% prctiel:              colour scale percentile: [0 - 100]
%
% OUTPUT:
% fig_knl:              figure with the kernels
%
% -- Nienke Blom, 31 March 2015

    %% PREPARATION
    input_parameters;
    
    [X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);
    set_figure_properties_bothmachines;

    [K2, parametrisation, which_wavefields, sameorownpercentile, prctiel] = checkargs(K, varargin(:));

    if strcmp(parametrisation, 'rhovsvp')
        params = {'rho2', 'vs2', 'vp2'};
    elseif strcmp(parametrisation, 'rhomulambda')
        params = {'rho', 'mu', 'lambda'};
    end

% disp(which_wavefields)
    switch which_wavefields
        case {'PSV', 'SH', 'total'}
            ncols = 1;
            veldje{1} = which_wavefields;
            % check if this field actually exists, otherwise error
%             params = fieldnames(K2);
            if ( ~isfield(K2.(params{1}), (veldje{1})) && ...
               ~isfield(K2.(params{2}), (veldje{1})) && ...
               ~isfield(K2.(params{3}), (veldje{1})) ); error('wavefield not existing'); 
            end
        case 'both'
%             ncols = 2;
%             veldje{1} = 'PSV';
%             veldje{2} = 'SH';
            % check if this field actually exists, otherwise error
%         case 'all'
            ncols = 3;
            veldje{1} = 'total';
            veldje{2} = 'PSV';
            veldje{3} = 'SH';
            % check if this field actually exists, otherwise error
%             params = fieldnames(K2);
%             if ( ~exist(K2.(params{1}).(veldje{1})) && ...
%                ~exist(K2.(params{2}).(veldje{1})) && ...
%                ~exist(K2.(params{3}).(veldje{1})) ); error('wavefield not existing'); end
        otherwise
            error('wavefields input not recognised');
    end

    switch sameorownpercentile % determines whether every kernel is plotted 
                               % on the same scale or not
        case 'same'
% disp 'same percentile'
            cmaxtype = 'fixed';
            K_together = [];
%             params = fieldnames(K2);
            for ii = 1:length(params)
                for jj = 1:length(veldje)
                    if isfield(K2.(params{ii}), veldje{jj})
                        K_together = [K_together; K2.(params{1}).(veldje{jj})];
%                         prctile(abs(K_together(:)),prctiel)
                    end
                end
            end
            scale = prctile(abs(K_together(:)),prctiel);
        case 'own'
            cmaxtype = 'perc';
            scale = prctiel;
    end
% disp(['cmaxtype: ',cmaxtype, ' - cmax: ', num2str(scale)])


    %% ACTUAL PLOTTING
    
    fig_knl = figure;
    set(fig_knl,'OuterPosition',pos_knl)
    % P-SV or 'total'
    for ii = 1:length(params); % loop over rows = parameters
        for jj = 1:ncols % loop over columns = which_waves fields ('PSV', 'SH', 'total')
            if isfield(K2.(params{ii}), veldje{jj})
                Ksm = filter_2Dfield(K2.(params{ii}).(veldje{jj}),smoothgwid);
                subplot(3,ncols,(ii-1)*ncols + jj)
                plot_kernel(X,Z,Ksm, ...
                    [params{ii}, ' - ', veldje{jj}],cmaxtype,scale);
            end
        end
    end
    
    % plotting the colorbar from -max to max
    h = colorbar('horiz');
    lmt = get(h, 'Limits');
    set(h,'XTick',[lmt(1) 0 lmt(2)])
    set(h,'Xticklabel',{'-max','0', 'max'})
    set(h,'Position',[0.3 0.05 0.4 0.02])
    
end

function [K2, parametrisation, which_waves, sameorownpercentile, prctiel] = checkargs(K, args)
    
    % determines the knl to be plotted, which waves (PSV, SH, both), the
    % prctile)

%     standard_prctile = 99.9;
    K2 = K;
    parametrisation = 'rhomulambda';
    which_waves = 'total';
    sameorownpercentile = 'same'; % 'same' or 'own'
    prctiel = 99.9;

    nargs = length(args);
    
    switch nargs
        case 0
% disp 'no args'
        case {1, 2, 3, 4, 5}
            Model = update_model();
            switch nargs
                case 1
% disp '1 arg'
                    if ischar(args{1}); 
                        which_waves = args{1};
                    else error('input to plot_kernels seems wrong');
                    end
                case 2
% disp '2 args'
                    if isstruct(args{2})
                        parametrisation = args{1};
                        Model = args{2};
                    elseif isnumeric(args{2})
                        sameorownpercentile = args{1};
                        prctiel = args{2};
                    else error('input to plot_kernels seems wrong');
                    end
                case 3
% disp '3 args'
                    if isstruct(args{2})
% disp 'a en b'
                        parametrisation = args{1};
                        Model = args{2};
                        which_waves = args{3};
                    elseif isnumeric(args{3})
% disp 'b en c'
                        which_waves = args{1};
                        sameorownpercentile = args{2};
                        prctiel = args{3};
                    else error('input to plot_kernels seems wrong');
                    end
                    
                case 4
% disp '4 args'
                    parametrisation = args{1};
                    Model = args{2};
                    sameorownpercentile = args{3};
                    prctiel = args{4};
                case 5
% disp '5 args'
                    parametrisation = args{1};
                    Model = args{2};
                    which_waves = args{3};
                    sameorownpercentile = args{4};
                    prctiel = args{5};
            end

            K2 = change_parametrisation_kernels('rhomulambda',parametrisation, K, Model);
            otherwise; error('input to plot_kernels seems wrong');
    end
    

end

