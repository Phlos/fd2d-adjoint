function [fig_waves, options, maks_normal, maks_diff] = video_waveprop(u_fw, varargin)

% show a video of wave propagation (and differential wave propagation)
% fig_waves = video_waveprop(u_fw, t_max)
%             will result in a single wave prop figure
% fig_waves = video_waveprop(u_fw, u_fw2, t_max, plotoptions)
%             will result in a figure with two subplots:
%             upper subplot: wave propagation u_fw
%             lower subplot: wave propabation u_fw - u_fw2



% initialise
input_parameters;
disp(['Assuming tmax as ',num2str(tmax),' seconds!']);
[X_new,Z_new,dx,dz]=define_computational_domain(Lx/1000,Lz/1000,nx,nz);
ntime = size(u_fw.x, 1);
[options, u_fw2, v_rec, v_rec2] = checkargs(varargin(:));

% options
% pause(1);

if options.make_movie
    warning('WARNING: making movie!')
end

% prepare figure
fig_waves = figure;
% pause(0.1);
if options.nplots == 2
    if options.plot_seismogram
        set(fig_waves, 'Position',[840,1138,769,1009])
        subplot(5,1,[1 2]);
        subplot(5,1,[3 4]);
        subplot(5,1,5);
    else
        set(fig_waves, 'Position', [850,1175,816,936]);
        subplot(2,1,1);
        subplot(2,1,2);
    end
else
    set(fig_waves, 'Position', [1100,1080,578,1085]);
end
pause(0.1);

% set some plotting options
nplots = options.nplots;
t_begin = options.t_begin;
t_end = options.t_end;
dt = tmax / ntime;
tee = dt : dt : tmax;

% calculate velocity seismogram
rec = options.plot_which_seis;
if options.plot_seismogram
    nufw = size(u_fw.z, 1);
    nvel = length(v_rec{rec}.z);
    sfe = nvel/nufw;
    vel = interp1(v_rec{rec}.z, [sfe : sfe : length(v_rec{rec}.z)]);
    if max(vel) < 0.1
        vel = vel*1000; unit = 'mm/s';
        if max(vel) < 0.1
            vel = vel*1000; unit = '{\mu}m/s';
        end
    end
    plot(vel); eilim = ylim;
    % also calculate vel2 if possible
    if ~(isnumeric(v_rec2) && isnan(v_rec2))
        vel2 = interp1(v_rec2{rec}.z, [sfe : sfe : length(v_rec{rec}.z)]);
        if max(vel2) < 0.1
            vel2 = vel2*1000; unit = 'mm/s';
            if max(vel2) < 0.1
                vel2 = vel2*1000; unit = '{\mu}m/s';
            end
        end
        hold on;
        plot(vel2); eilim2 = ylim;
        eilim(1) = min(eilim(1), eilim2(1));
        eilim(2) = max(eilim(2), eilim2(2));
    end
end

% loop to obtain wave propagation from supplied info
for itime = t_begin:t_end
    
    %% calculation / preparation
    % time
    time = itime/ntime * tmax;
    
    %- calculate 'standard' wave prop plot
    ux = squeeze(u_fw.x(ntime - itime + 1,:,:));
    uz = squeeze(u_fw.z(ntime - itime + 1,:,:));
    u_mag = sqrt(ux.^2 + uz.^2);
    maks_normal(itime) = max(abs(u_mag(:)));
    
    %     % divergence of u
    %     [duxdx, ~] = gradient(ux, dx);
    %     [~, duzdz] = gradient(uz, dz);
    %     udiv = duxdx;% + duzdz;
    
    % calculate differential wave prop plot
    if options.plot_diff_wave_prop
        ux_diff = squeeze(u_fw.x(ntime - itime + 1,:,:)) - ...
            squeeze(u_fw2.x(ntime - itime + 1,:,:));
        uz_diff = squeeze(u_fw.z(ntime - itime + 1,:,:)) - ...
            squeeze(u_fw2.z(ntime - itime + 1,:,:));
        u_mag_diff = sqrt(ux_diff .^2 + uz_diff .^2);
        maks_diff(itime) = max(abs(u_mag_diff(:)));
    end
    
    %% actual plotting
    clf
    
    %- plot 'standard'
    if options.plot_seismogram
        subplot1 = subplot(5,1, [1 2]);
    else
        subplot1 = subplot(2, 1, 1);
    end
    pcolor(X_new,Z_new, u_mag');
    pl_opt = options;
    if isfield(options, 'movie_label1')
        pl_opt.movie_label = options.movie_label1;
    end
    make_subplot_shiny(subplot1, time, pl_opt);
    caxis_normal = caxis;
    
    %- plot differential
    if options.plot_diff_wave_prop
        if options.plot_seismogram
            subplot2 = subplot(5,1, [3 4]);
        else
            subplot2 = subplot(nplots, 1, 2);
        end
        pcolor(X_new,Z_new, u_mag_diff');
        pl_opt = options;
        pct = pl_opt.pct;
        pl_opt.max_text    = ['amplitude ', num2str(pct*100, '%1.0f'), '%'];
        pl_opt.color_lim   = caxis_normal*pct;
        if isfield(options, 'movie_label2')
            pl_opt.movie_label = options.movie_label2;
        end
        make_subplot_shiny(subplot2, time, pl_opt);
    end
    
    %- plot seismogram
    if options.plot_seismogram
        subplot3 = subplot(5,1,5);
        hold on;
        % plot the reference wave prop, if existant
        if exist('vel2', 'var')
            plot(tee(1:itime), vel2(1:itime), '-r');
        end
        % plot the actual (with anomaly) seismogram.
        plot(tee(1:itime), vel(1:itime), '-k'); 
        plot(tee(itime), vel(itime), 'ro');
        xlim([0 tmax]);
%         maks_ei = max(abs(vel)); ylim([-maks_ei maks_ei]);
        ylim(eilim);
        xlabel('time [s]'); ylabel(['v (rec) [', num2str(unit),']'])
        
        legend('homogeneous', 'model with anomaly', ...
            'Location', 'SouthWest');
        
    end
    
    %% get movie frame
    if options.make_movie
        if exist('movie_index','var')
            movie_index=movie_index+1;
        else
            pause(10);
            movie_index=1;
        end
        M(movie_index)=getframe(gcf);
    end
    
    pause(0.01)
    
end


%% movie
if options.make_movie
    save('./output/movie_M.mat', 'M');
    clearvars('u_fw*', 'varargin');
    make_the_movie(M, options);
end

end

function [options, u_fw2, v_rec, v_rec2] = checkargs(args)

% set plotting options for video_waveprop

% initialise
input_parameters;

% standard values
options.plot_diff_wave_prop = false;
options.plot_seismogram = false;
options.pct = 0.20;
options.nplots = 1;
options.t_begin = 1;
options.t_end = tmax;
options.make_movie = false;
options.movie_file = './output/differential_wave_prop';
options.wavefield_type = 'displacement';
u_fw2 = NaN;
v_rec = NaN;
v_rec2 = NaN;


for iarg = 1:numel(args)
    if isstruct(args{iarg})
        if(isfield(args{iarg}, 'x'))
            u_fw2 = args{iarg};
            options.plot_diff_wave_prop = true;
            options.nplots = options.nplots +1;
        else
            % loop over the fields in the 'plotoptions' struct.
            fn = fieldnames(args{iarg});
            for ii = 1:numel(fn)
                options.(fn{ii}) = args{iarg}.(fn{ii});
            end
        end
        
%     elseif isnumeric(args{iarg})
%         
%         % the argument is not part of a structure and numeric, so t_max
%         t_max = args{iarg};
        
    elseif iscell(args{iarg})
        % we've got v_rec here now
        if iscell(v_rec);
            v_rec2 = args{iarg};
        else
            v_rec = args{iarg};
        end
        
        
    else
        error('argument is not a struct');
    end
    
end

if isnumeric(v_rec) && isnan(v_rec)
    warning('not plotting a seismogram because we don''t have v_rec!')
    options.plot_seismogram = false;
% else
%     nufw = size(u_fw.x, 1);
%     nvel = size(v_rec{1}.x);
%     sfe = nvel/nufw;
%     vel = interp1(v_rec{1}.x, [500);
end

if options.t_begin < 1
    options.t_begin = 1;
end

if (options.t_end > tmax)
    options.t_end = tmax;
end

if isfield(options, 'Model_rho') || isfield(options, 'Model_norho')
    if ~isfield(options, 'Model_rho')
        options.Model_rho.rho = zeros(nx,nz);
    end
    if ~isfield(options, 'Model_norho')
        options.Model_rho.rho = zeros(nx,nz);
    end
    options.mod_rho_diff = options.Model_rho.rho - options.Model_norho.rho;
end

if ~isfield(options, 'plot_which_seis')
    options.plot_which_seis = 1;
    disp 'WARNING! assuming first receiver''s seismogram is the one you want'
end

end

function ssplot= make_subplot_shiny(ssplot, time, options)

    %% initialise
    input_parameters;
    [X_new,Z_new,~,~] = define_computational_domain(Lx/1000, Lz/1000, nx, nz);

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
    
%     colorbar;
    
    %% shading, axis
    shading flat; axis image;
    hold on;
    xlabel('x [km]');
    ylabel('z [km]');
    
    %% plot src & rec
    for k=1:length(src_info)
        src_x_new(k) = src_info(k).loc_x / 1000;
        src_z_new(k) = src_info(k).loc_z / 1000;
    end
    plot(src_x_new,src_z_new,'kx')
    
    rec_x_new = rec_x ./ 1000;
    rec_z_new = rec_z ./ 1000;
    plot(rec_x_new,rec_z_new,'ko')
    
    %% plot contour
    if isfield(options, 'mod_rho_diff')
        contour(X_new,Z_new, options.mod_rho_diff', 1, '--k');
    end
    
    %% plot timestamp
    timestamp_text=['t = ',num2str(floor(time),'%4.0f'), ' s'];
    timestamp_plot = text(0.05,0.92, timestamp_text, ...
        'HorizontalAlignment','left', 'Units', 'normalized');
    
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
    text(0.95,0.92, max_text, ...
        'HorizontalAlignment','right', 'Units', 'normalized');
    
    %% plot copyright
    text(0.05,0.08, ['\copyright Nienke Blom'], 'Color', [0.72 0.72 0.72], ...
        'HorizontalAlignment','left', 'Units', 'normalized');
    
    %% plot wavefield description
    if isfield(options, 'movie_label')
        movie_label = options.movie_label;
    end
    if exist('movie_label', 'var')
    text(0.95,0.08, movie_label, ...
        'HorizontalAlignment','right', 'Units', 'normalized');
    end

end

function make_the_movie(M, options)

    movie_file = options.movie_file;
    disp(['storing movie: ',movie_file, '.avi']);
    writerObj=VideoWriter([movie_file, '.avi'],'Motion JPEG AVI');  
    writerObj.FrameRate = 10;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
    % frame rate 20 fps:
    writerObj=VideoWriter([movie_file,'.20-fps.avi'],'Motion JPEG AVI');
    writerObj.FrameRate = 20;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);

end

% function seis = plot_seism(u_fw, u_fw2, itime)
% 
%     % plot seismogram in video
%     
%     % initialise
%     input_parameters;
%     [~,~,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
%     [rec_x_id,rec_z_id,n_receivers] = compute_indices(rec_x,rec_z,Lx,Lz,dx,dz);
%     
%     seis(itime)
%     
%     % plot
%     
% 
% end
