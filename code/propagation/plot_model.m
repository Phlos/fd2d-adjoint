function fig_mod = plot_model(varargin)


% function that plots the model in rho mu lambda parametrisation.
% NOTE: Currenly this is a rather ugly function and I could make it much
% nicer by putting the plotting commands within a for loop over rho mu
% lambda.
%
% INPUT:
% Model:    struct containing .rho .mu .lambda
%
% OUTPUT:   
% - figure with the model plotted. Colour scale are the actual max and min
%   of the parameter values, not some standard deviation.

% format long

[Model, middle] = checkargs(varargin);


input_parameters;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
set_figure_properties_maggi;
% mu = Model.mu;
% rho = Model.rho;
% lambda = Model.lambda;

load 'propagation/cm_model.mat';

fig_mod = figure;
set(fig_mod,'OuterPosition',pos_mod) 
% set(gca,'FontSize',14)

j=1;
for params = fieldnames(Model)';
    

    param = Model.(params{1});

    subplot(1,3,j)
    
    hold on
    pcolor(X,Z,param');
    
    % difference_mumax_mumin = max(mu(:)) - min(mu(:));
    
    %- colour scale
    if (middle == 1234567890)
        if all(param == param(1))
            cmax = param(1) + 0.01*param(1);
            cmin = param(1) - 0.01*param(1);
            %     disp 'bips!!! all param are the same!'
        else
            % max and min are calculated in this way so that the most common value
            % (i.e. the background value) is white, and that the extreme coulours
            % are determined by whichever of max and min is farthest off from the
            % background colour.
            %     cmid = mode(param(:))
            [bincounts, centre] = hist(param(:),100);
            [~,ix]=max(bincounts);
            cmid = centre(ix);
            cmax = cmid + max(abs(param(:) - cmid));
            cmin = cmid - max(abs(param(:) - cmid));
            %     cmax = max(param(:));
            %     cmin = 2*cmid - cmax;
            %     disp 'the param are not all the same'
        end
    else
            % max and min are calculated in this way so that the most common value
            % (i.e. the background value) is white, and that the extreme coulours
            % are determined by whichever of max and min is farthest off from the
            % background colour.
            %     cmid = mode(param(:))
%             [bincounts, centre] = hist(param(:),100);
%             [~,ix]=max(bincounts);
            cmid = middle;
            cmax = cmid + max(abs(param(:) - cmid));
            cmin = cmid - max(abs(param(:) - cmid));
            %     cmax = max(param(:));
            %     cmin = 2*cmid - cmax;
            %     disp 'the param are not all the same'
    end
    caxis([cmin cmax]);
    
    
    for k=1:length(src_x)
        plot(src_x(k),src_z(k),'kx')
    end
    
    for k=1:length(rec_x)
        plot(rec_x(k),rec_z(k),'ko')
    end
    colormap(cm_model);
    axis image
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
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    hold off;
    
    j = j+1;
end

end

function [Model, middle, outparam] = checkargs(arg)

% checks the input argument and defines the fields to be plotted, the
% middle of the plot colour (in all three) and 

narg = length(arg);

Model = arg{1};

switch narg
    case 3
%         disp '3 arguments'
        middle = arg{2};
        outparam = arg{3};
        Model = change_parametrisation('rhomulambda','rhovsvp',Model);
    case 2
%          disp '2 arguments'

        if ischar(arg{2})
            
            outparam = arg{2};
            Model = change_parametrisation('rhomulambda','rhovsvp',Model);
            
            middle = 1234567890;
            
        elseif isnumeric(arg{2})
            middle = arg{2};
            teststruct.mu = 1;
            teststruct.lambda = 2;
            teststruct.rho = 3;
            Model = orderfields(Model, teststruct);
        end
        

    case 1
%         disp '1 argument'
        middle = 1234567890;
    otherwise
        error('wrong number of variable arguments')
end

end