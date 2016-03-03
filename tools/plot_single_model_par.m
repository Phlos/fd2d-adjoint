function fig_par = plot_single_model_par(Model, Model_ref, which_param, plotoptions)
    % plots a single model parameter
    
    plotoptions = check_options(plotoptions);
    
    %% preparation
    input_parameters;
    fig_par = figure;
    set(fig_par, 'OuterPosition', [1 413 1059 668]);
    
    load 'code/propagation/cm_model.mat';
    
    %% recalculate lengths
    
    [X,Z,~,~]=define_computational_domain(Lx,Lz,nx,nz);
    [X,Z, srcs, recs] = recalculate_to_km(X, Z, src_info, rec_x, rec_z);

    %% determine parameter to be plotted
    
    if any(strcmp(which_param, {'vs', 'vp'}))
        Model = change_parametrisation('rhomulambda', 'rhovsvp', Model);
        Model_ref = change_parametrisation('rhomulambda', 'rhovsvp', Model_ref);
    end
    
    param = Model.(which_param) - Model_ref.(which_param);
    
    %% plotting
    
    % make sure that the Y axis points downwards
    set(gca, 'YDir', 'reverse');
    
    hold on
    colorplot = pcolor(X,Z,param');
    
    % difference_mumax_mumin = max(mu(:)) - min(mu(:));
    
    if plotoptions.cscale_own
    %- colour scale
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
        caxis([cmin cmax]);
    else
        caxis(plotoptions.cscale_minmax.(which_param));
    end
    
    colormap(cm_model);
    axis image
    shading flat
    
    % plot sources and receivers
    if plotoptions.plotsrcrec
        plot(srcs.x,srcs.z,'kx','LineWidth',2,'MarkerSize',8)
        plot(recs.x,recs.z,'ko','LineWidth',2,'MarkerSize',8)
    end
    
    % plot contours of rho-vs-vp anomalies

    if plotoptions.plotrvvcontours
        Mod_rvv = change_parametrisation('rhomulambda', 'rhovsvp', Model);
        %     MR = update_model(85);
        MR = plotoptions.Mod_contours;
        MBG = plotoptions.Mod_contours_bg;
        Model_real_rvv = change_parametrisation('rhomulambda', 'rhovsvp', MR);
        Mod_ref_rvv = change_parametrisation('rhomulambda', 'rhovsvp', MBG);
        %         if ~strcmp(which_param, 'rho')
            contour(X, Z, Model_real_rvv.rho' - Mod_ref_rvv.rho', 2, '--k');
%         end
%         if ~strcmp(which_param, 'vs')
            contour(X, Z, Mod_rvv.vs' - Mod_ref_rvv.vs', 2, '--k');
%         end
%         if ~strcmp(which_param, 'vp')
            contour(X, Z, Mod_rvv.vp' - Mod_ref_rvv.vp', 2, '--k');
%         end
    end
    

    if plotoptions.plotaxes
        switch which_param
            case 'rho'
                title([which_param, ' [kg/m^3]'])
            case {'vs', 'vp'}
                title([which_param, ' [m/s]'])
            case {'mu', 'lambda'}
                title([which_param, ' [N/m^2]'])
            otherwise
                title([which_param, ' [unit??]'])
        end
        xlabel('x [km]');
        ylabel('z [km]');
        cbar = colorbar;
        %     ylabel(cbar, '
    end
    hold off;
    
    if ~plotoptions.plotaxes
        axis off
    end
    
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

function plotoptions = check_options(plotoptions)
    
    input_parameters;
    
    if ~isfield(plotoptions, 'plotaxes')
        plotoptions.plotaxes = true;
    end
    
    if ~isfield(plotoptions, 'plotsrcrec')
        plotoptions.plotsrcrec = true;
    end
    
    if ~isfield(plotoptions, 'plotrvvcontours')
        plotoptions.plotrvvcontours = true;
    end
    
    if ~isfield(plotoptions, 'Mod_contours')
        plotoptions.Mod_contours = update_model(true_model_type);
    end
    
    if ~isfield(plotoptions, 'Mod_contours_bg')
        plotoptions.Mod_contours_bg = update_model(bg_model_type);
    end
    
    if ~isfield(plotoptions, 'cscale_own')
        plotoptions.cscale_own = false;
    end
    
    if ~isfield(plotoptions, 'cscale_minmax')
        mmrho = 53.2870;
        mmvs  = 71.1881;
        mmvp  = 132.9489;
        mmmu      = 1.8932e10;
        mmlambda  = 5.4279e9;
        plotoptions.cscale_minmax.rho = [-mmrho mmrho];
        plotoptions.cscale_minmax.vs  = [-mmvs mmvs];
        plotoptions.cscale_minmax.vp  = [-mmvp mmvp];
        plotoptions.cscale_minmax.mu  = [-mmmu mmmu];
        plotoptions.cscale_minmax.lambda  = [-mmlambda mmlambda];
    end
    
end