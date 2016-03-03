function fig_mod = plot_model(varargin)


% function that plots the model in rho mu lambda parametrisation.
% NOTE: Currenly this is a rather ugly function and I could make it much
% nicer by putting the plotting commands within a for loop over rho mu
% lambda.
%
% - fig_mod = plot_model(Model);
% - fig_mod = plot_model(Model, outparam);
% - fig_mod = plot_model(Model, middle);
% - fig_mod = plot_model(Model, middle, outparam);
%
% INPUT:
% Model:    struct containing .rho .mu .lambda -OR- modelnr.
% outparam: string, 'rhomulambda', 'rhovsvp' (future: 'rhomukappa'?)
%           parametrisation of output plot (also parametrisation of middle)
% middle:   1x3 array w/ colour scale centre values for [param1, param2, param3]; 
%
% OUTPUT:   
% - figure with the model plotted. Colour scale are defined by middle (if 
%   given) -- otherwise, it is the actual max and min of the parameter 
%   values, not some standard deviation.


[Model, middle] = checkargs(varargin);

if ~exist('plot_UM_separate', 'var')
    plot_UM_separate = true;
end

%% preparation
input_parameters;

set_figure_properties_bothmachines;
if plot_UM_separate
    pos_mod = pos_mod .* [1 0.75 1.5 1.5];
end

load 'propagation/cm_model.mat';


%% recalculation of lengths
[X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);

[X,Z, srcs, recs] = recalculate_to_km(X, Z, src_info, rec_x, rec_z);


%% figure
% fig_mod = figure('Visible', 'off');
fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod)
if (feature('showfigurewindows') == 0)
    set(fig_mod, 'PaperUnits', 'points');
    set(fig_mod, 'PaperPosition', pos_mod);
end

j=1;
for params = fieldnames(Model)';
    

    param = Model.(params{1});
%     max(abs(param(:)))

%     g = subplot(2,3,j);
%     p = get(g,'position');
% %     p(4) = p(4)*1.50;  % Add 10 percent to height
%     set(g, 'position', p);

    %% plot whole mantle
    
    if plot_UM_separate
        subplot(2,3,3+j);
    else
        subplot(1,3,j);
    end
%     subplot(4,3,[3+j 6+j]);

    % make sure that the Y axis points downwards
    set(gca, 'YDir', 'reverse');
    
    hold on
    pcolor(X,Z,param');
    
    % difference_mumax_mumin = max(mu(:)) - min(mu(:));
    
    %- colour scale
%     disp 'bips'
    if (isnan( middle.(params{1}) ))
%         param_min_mode_param = max(abs(param(:)-mode(param(:))))
%         mode_param_1 = mode(param(1))
        if all(abs(param-mode(param(:))) <= max(1e-10*mode(param(1)), 1e-10))
            cmax = param(1) + max(0.01*param(1), 1);
            cmin = param(1) - max(0.01*param(1), 1);
%                 disp 'bips!!! all param are the same!'
        else
            % max and min are calculated in this way so that the most common value
            % (i.e. the background value) is white, and that the extreme coulours
            % are determined by whichever of max and min is farthest off from the
            % background colour.
            % this is done not with mode but with hist, because after
            % update, the values all vary just a bit.
            %     cmid = mode(param(:))
            [bincounts, centre] = hist(param(:),100);
            [~,ix]=max(bincounts);
            cmid = centre(ix);
            cmax = cmid + max(abs(param(:) - cmid));
            cmin = cmid - max(abs(param(:) - cmid));

        end
        
    else
        
        refval = param(1);
        if refval == 0; refval = 1; end
        if all(abs(param-mode(param(:))) <= 1e-10*mode(param(:)))
            cmax = param(1) + 0.01*param(1);
            cmin = param(1) - 0.01*param(1);
%                 disp 'bips!!! all param are the same!'
        else
            
            % max and min are calculated in this way so that the the 
            % background value) is white, and that the extreme coulours are
            % determined by whichever of max and min is farthest off from 
            % the background colour.
            
            cmid = middle.(params{1});
            cmax = cmid + max(abs(param(:) - cmid));
            cmin = cmid - max(abs(param(:) - cmid));
            
        end
        
    end
    
    caxis([cmin cmax]);
    
    % plot sources and receivers
    plot(srcs.x,srcs.z,'kx','LineWidth',0.3,'MarkerSize',4)
    plot(recs.x,recs.z,'ko','LineWidth',0.3,'MarkerSize',4)

    colormap(cm_model);
    axis image
    shading flat
    if ~plot_UM_separate
        switch params{1}
            case 'rho'
                title([params{1}, ' [kg/m^3]'])
            case {'vs', 'vp'}
                title([params{1}, ' [m/s]'])
            case {'mu', 'lambda'}
                title([params{1}, ' [N/m^2]'])
            otherwise
                title([params{1}, ' [unit??]'])
        end
    end
    xlabel('x [km]');
    ylabel('z [km]');
    colorbar
    hold off;
    
    % plot description
    if plot_UM_separate
        text(0.5, 0.1, ['whole mantle'], ...
            'Units', 'normalized', 'HorizontalAlignment','center');
    end
    
    %% plot UM zoom in
%     subplot(4,3,j)
if plot_UM_separate
    subplot(2,3,j);
    set(gca, 'YDir', 'reverse');
    hold on
    pcolor(X,Z,param');
    caxis([cmin cmax]);
    % plot sources and receivers
    plot(srcs.x,srcs.z,'kx','LineWidth',0.3,'MarkerSize',4)
    plot(recs.x,recs.z,'ko','LineWidth',0.3,'MarkerSize',4)
    
    colormap(cm_model);
    %     axis image
    shading flat
    switch params{1}
        case 'rho'
            title([params{1}, ' [kg/m^3]'])
        case {'vs', 'vp'}
            title([params{1}, ' [m/s]'])
        case {'mu', 'lambda'}
            title([params{1}, ' [N/m^2]'])
        otherwise
            title([params{1}, ' [unit??]'])
    end
    ylim([0 660]);
    xlim([min(X(:)), max(X(:))]);
    %     xlabel('x [km]');
    ylabel('z [km]');
    colorbar
    
    % plot description
    text(0.5, 0.1, ['upper mantle (vertical stretch)'], ...
        'Units', 'normalized', 'HorizontalAlignment','center');
end
    
%     % plot histogram
%     subplot(4,3,6+j)
%     
%      hist(Model.(params{1})(:),100)
%      h = findobj(gca,'Type','patch');
%      set(h,'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    

    %% next plot
    j = j+1;
    
%     disp(['cmin: ',num2str(cmin,'%5.5e'),'   cmax: ', num2str(cmax,'%5.5e')]);
%     cmax
end

end

function [Model, middle, outparam] = checkargs(arg)

% checks the input argument and defines the fields to be plotted, the
% middle of the plot colour (in all three) and 

narg = length(arg);


% % loop over input arguments
% for ii = 1:narg
%     if isstruct(arg{ii})
%         Model = arg{ii};
%         if length(Model) > 1
%             error('You supplied more than one models (i.e. a struct w/ multiple models?)');
%         end
%     elseif isnumeric(arg{ii}) && numel(arg{ii}) == 1
%         modelnr = arg{ii};
%         Model = update_model(modelnr);
%     elseif isnumeric(arg{ii}) && numel(arg{ii}) == 3
%         middle = arg{ii};
%     elseif islogical(arg{ii})
%         plot_UM_separate = arg{ii};
%     end
% end

switch narg
    
    case 1
        % disp '1 argument'
        
        if isstruct(arg{1})
            Model = arg{1};
            if length(Model) > 1
                error('You supplied more than one models (i.e. a struct w/ multiple models?)');
            end
            
        elseif isnumeric(arg{1})
            modelnr = arg{1};
            Model = update_model(modelnr);
        else
            error('input has to be a model structure or modelnr');
        end
        %         middle = [NaN NaN NaN];
        params = fieldnames(Model);
        for ii = 1:length(params)
            middle.(params{ii}) = NaN;
        end
%         middle.rho = NaN;
%         middle.mu = NaN;
%         middle.lambda = NaN;
        
    case 2
        % disp '2 arguments'
        if isstruct(arg{1})
        Model = arg{1};
            if length(Model) > 1
                error('You supplied more than one models (i.e. a struct w/ multiple models?)');
            end
        elseif isnumeric(arg{1})
            modelnr = arg{1};
            Model = update_model(modelnr);
        else
            error('1st input has to be a model or a modelnr');
        end
        

        if ischar(arg{2})
            
            outparam = arg{2};
            Model = change_parametrisation('rhomulambda',outparam,Model);
            
            % middle = [NaN NaN NaN];
            middle1.rho = NaN;
            middle1.mu = NaN;
            middle1.lambda = NaN;
            middle = change_parametrisation('rhomulambda',outparam,middle1);
            
        elseif isstruct(arg{2})
            middle = arg{2};
        else
            disp 'allowed var types for input argument {2}: char, struct'
            error('the var type of input argument {2} was not recognised')
        end
        
    case 3
        % disp '3 arguments'
        Model = arg{1};
        middle1 = arg{2};
        outparam = arg{3};
        
        if length(Model) > 1
            error('You supplied more than one models (i.e. a struct w/ multiple models?)');
        end
        
        Model = change_parametrisation('rhomulambda',outparam,Model);
        if middle1.rho == 0 && middle1.mu == 0 && middle1.lambda == 0
            if strcmp(outparam, 'rhovsvp')
                middle.rho = 0; middle.vs = 0; middle.vp = 0;
            elseif strcmp(outparam, 'rhomulambda')
                middle.rho = 0; middle.mu = 0; middle.lambda = 0;
            else
                error('unknown output parametrisation');
            end
        else
            middle = change_parametrisation('rhomulambda',outparam,middle1);
        end

    otherwise
        error('wrong number of variable arguments')
end

% fn = fieldnames(middle);
% for ii = 1:3
%     if isnan(middle.(fn{ii}))
%         middle.(fn{ii}) = 0;
%     end
% end

end

function [X, Z, srcs, recs] = recalculate_to_km(X, Z, src_info, rec_x, rec_z);
    
    % convert distances to km & set Z to depth below surface
surface_level = 2890; % only valid in PREM
X = X ./ 1000;
Z = Z ./ 1000;
Z = surface_level - (Z);
for k=1:length(src_info)
    src_x(k) = src_info(k).loc_x / 1000;
    src_z(k) = src_info(k).loc_z / 1000;
    src_z(k) = surface_level - src_z(k);
end

for k=1:length(rec_x)
    rec_x(k) = rec_x(k) ./ 1000;
    rec_z(k) = rec_z(k) ./ 1000;
    rec_z(k) = surface_level - rec_z(k);
end

srcs.x = src_x;
srcs.z = src_z;

recs.x = rec_x;
recs.z = rec_z;

end