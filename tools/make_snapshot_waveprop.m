function [fig_snapshot, options] = make_snapshot_waveprop(time, u_fw, varargin)

% make snapshots of wave propagation file
% fig_waves = make_snapshot_waveprop(time, u_fw, t_max)
%             will result in a single wave prop figure
% fig_waves = make_snapshot_waveprop(time, u_fw, u_fw2, t_max, plotoptions)
%             will result in a figure with two subplots:
%             upper subplot: wave propagation u_fw
%             lower subplot: wave propabation u_fw - u_fw2



% initialise
input_parameters;
[X_new,Z_new,dx,dz]=define_computational_domain(Lx/1000,Lz/1000,nx,nz);
ntime = size(u_fw.x, 1);
[options, u_fw2, v_rec, fig_snapshot] = checkargs(varargin(:));

% options
% pause(1);

% prepare figure
if ~isgraphics(fig_snapshot)
    fig_snapshot = figure;
end
% pause(0.1);

    if options.plot_seismogram
        set(fig_snapshot, 'Position', [1100,1080,578,1085]);
    else
        set(fig_snapshot, 'Position', [426,1352,1248,712]);
    end

pause(0.1);

% set some plotting options
dt = tmax / ntime;
tee = dt : dt : tmax;



%% calculation / preparation
% time
itime = round(time / tmax * ntime);

% calculate normal wave prop
ux = squeeze(u_fw.x(ntime - itime + 1,:,:));
uz = squeeze(u_fw.z(ntime - itime + 1,:,:));
u_mag = sqrt(ux.^2 + uz.^2);
maks_normal(itime) = max(abs(u_mag(:)));

% calculate differential wave prop
if options.plot_diff_wave_prop
    ux_diff = squeeze(u_fw.x(ntime - itime + 1,:,:)) - ...
        squeeze(u_fw2.x(ntime - itime + 1,:,:));
    uz_diff = squeeze(u_fw.z(ntime - itime + 1,:,:)) - ...
        squeeze(u_fw2.z(ntime - itime + 1,:,:));
    u_mag_diff = sqrt(ux_diff .^2 + uz_diff .^2);
    maks_diff(itime) = max(abs(u_mag_diff(:)));
else

end

%% actual plotting
clf

if options.plot_seismogram
    subplot1 = subplot(3,1, [1 2]);
else
    subplot1 = subplot(1,1,1);
end

%- plot differential
if options.plot_diff_wave_prop

%     pcolor(X_new,Z_new, u_mag_diff');
    pl_opt = options;
    if isfield(options, 'pct')
        pct = options.pct;
    else
        pct = 0.20;
    end
    pl_opt.color_lim = [-maks_normal(itime) maks_normal(itime)]*pct;
    pl_opt.max_text    = ['amplitude ', num2str(pct*100, '%1.0f'), '%'];
    if isfield(options, 'text_label2')
        pl_opt.text_label = options.text_label2;
    end
    plot_wavefield_shiny(subplot1, u_mag_diff, time, pl_opt);
else
    pl_opt = options;
    if isfield(options, 'text_label1')
        pl_opt.text_label = options.text_label1;
    end
    plot_wavefield_shiny(subplot1, u_mag, time, pl_opt);
    caxis_normal = caxis;
end

%- plot seismogram
if options.plot_seismogram
    subplot3 = subplot(3,1,3);
    
    nufw = size(u_fw.z, 1);
    nvel = length(v_rec{1}.z);
    sfe = nvel/nufw;
    vel = interp1(v_rec{1}.z, [sfe : sfe : length(v_rec{1}.z)]);
    if max(vel) < 0.1
        vel = vel*1000; unit = 'mm/s';
        if max(vel) < 0.1
            vel = vel*1000; unit = '{\mu}m/s';
        end
    end
    plot(vel); eilim = ylim;
    
    plot(tee(1:itime), vel(1:itime), '-k'); hold on;
    plot(tee(itime), vel(itime), 'ro');
    xlim([0 tmax]);
    ylim(eilim);
    xlabel('time [s]'); ylabel(['v (rec) [', num2str(unit),']'])
    
end


end

function [options, u_fw2, v_rec, fig_snapshot] = checkargs(args)

% set plotting options for video_waveprop

% initialise
input_parameters;

% default values
options.plot_diff_wave_prop = false;
options.plot_seismogram = false;
options.plot_copyright = false;
options.plot_label = false;
options.plot_timestamp = true;
options.wavefield_type = 'displacement';
options.zoom_box = [0 0 Lx/1000 Lz/1000];
options.no_text = false;

% options.text_label1 = '';
% options.text_label2 = '';
u_fw2 = NaN;
v_rec = NaN;
fig_snapshot = NaN;


for iarg = 1:numel(args)
    if isstruct(args{iarg})
        if(isfield(args{iarg}, 'x'))
            u_fw2 = args{iarg};
            options.plot_diff_wave_prop = true;
        else
            fn = fieldnames(args{iarg});
            for ii = 1:numel(fn)
                options.(fn{ii}) = args{iarg}.(fn{ii});
            end
        end
        
%     elseif isnumeric(args{iarg})
%         tmax = args{iarg};
        
    elseif iscell(args{iarg})
        % we've got v_rec here now
        v_rec = args{iarg};
        
    elseif isgraphics(args{iarg})
        fig_snapshot = args{iarg};
        
    else
        error('argument is not a struct');
    end
    
end

if isnumeric(v_rec) && isnan(v_rec)
    options.plot_seismogram = false;
end


if isfield(options, 'Model_rho') || isfield(options, 'Model_norho')
    if ~isfield(options, 'Model_rho')
        options.Model_rho.rho = zeros(nx,nz);
    end
    if ~isfield(options, 'Model_norho')
        options.Model_norho.rho = zeros(nx,nz);
    end
    options.mod_rho_diff = options.Model_rho.rho - options.Model_norho.rho;
end

end

function ssplot= plot_wavefield_shiny(ssplot, u_mag, time, options)

    %% initialise
    input_parameters;
    [X_new,Z_new,~,~] = define_computational_domain(Lx/1000, Lz/1000, nx, nz);
    zoom_box = options.zoom_box;
    
    %% plot
    X_new = X_new - zoom_box(1);
    Z_new = Z_new - zoom_box(2);
    pcolor(X_new,Z_new, u_mag');
    
    

    %% colour scale
    load('code/propagation/cm_coolwarm.mat')
    colormap(colormap_CoolWarm)
%     load('code/propagation/cm_bluewhitered_soft.mat')
%     colormap(colormap_bluewhitered_soft)
    if (isfield(options, 'color_lim'))
        maks = options.color_lim(2);
        caxis(options.color_lim);
    else
        maks = max(abs(caxis));
        caxis([-maks maks]);
    end
    
    
    %% shading, axis
    shading interp; 
    axis image;
    hold on;
    xlim([0 zoom_box(3)-zoom_box(1)]);
    ylim([0 zoom_box(4)-zoom_box(2)]);
    xlabel('x [km]');
    ylabel('z [km]');
    
    %% plot src & rec
    for k=1:length(src_info)
        src_x_new(k) = src_info(k).loc_x / 1000 - zoom_box(1);
        src_z_new(k) = src_info(k).loc_z / 1000 - zoom_box(2);
    end
    plot(src_x_new,src_z_new,'kx', 'MarkerSize', 15);%, 'LineWidth', 3)
    
    rec_x_new = rec_x ./ 1000 - zoom_box(1);
    rec_z_new = rec_z ./ 1000 - zoom_box(2);
    plot(rec_x_new,rec_z_new,'ko', 'MarkerSize', 15);%, 'LineWidth', 3)
    
    %% plot contour
    if isfield(options, 'mod_rho_diff')
        contour(X_new,Z_new, options.mod_rho_diff', 1, '--k');
    end
    
    %% plot timestamp
    if options.plot_timestamp && ~options.no_text
        timestamp_text=['t = ',num2str(floor(time),'%4.0f'), ' s'];
        timestamp_plot = text(0.05,0.92, timestamp_text, ...
            'HorizontalAlignment','left', 'Units', 'normalized');
    end
    
    %% plot maxtext
    unit = ' m'; 
    if maks < 0.1
        maks = maks * 1000; unit = ' mm';
        if maks < 0.1
            maks = maks * 1000; unit = ' {\mu}m'; 
        end
    end
    if ( isfield(options, 'wavefield_type') && strcmp(options.wavefield_type, 'velocity') )
        unit = [unit,'/s']; 
    end
    if ( isfield(options, 'max_text') )
        max_text = options.max_text;
    else
        max_text = ['amplitude = \pm', num2str(maks,'%3.1f'), unit];
    end
    if  ~options.no_text
        text(0.95,0.92, max_text, ...
            'HorizontalAlignment','right', 'Units', 'normalized');
    end
    %% plot copyright
    if options.plot_copyright && ~options.no_text
        text(0.05,0.08, '\copyright Nienke Blom', 'Color', [0.72 0.72 0.72], ...
            'HorizontalAlignment','left', 'Units', 'normalized');
    end
    
    %% plot wavefield description
    if options.plot_label && ~options.no_text
        if isfield(options, 'text_label')
            text_label = options.text_label;
        end
        if exist('text_label', 'var')
            text(0.95,0.08, text_label, ...
                'HorizontalAlignment','right', 'Units', 'normalized');
        end
    end

end