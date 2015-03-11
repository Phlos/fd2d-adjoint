% check compatibility of model and modelling values.
%
% - the maximum frequency of the source etc.
% - whether dx and dz are the same
% - whether the absorbing boundaries aren't too wide for the domain size
% - stability condition
% - nt long enough for signal to reach receiver

vp = sqrt( (lambda + 2*mu) ./ rho );
vs = sqrt( mu ./ rho );

%%  stability condition
dt_max = 1/sqrt(2) * min(dx,dz)/ max(vp(:));
disp(['dt is ', num2str(dt), ' and dtmax is ', num2str(dt_max)]);

if (dt_max < dt)
    error({'ERROR: stability condition not satisfied! dt is too large'; ...
           ['dt needs to be at most ',num2str(dt_max)] } )
end

%% dx and dz same size -- 1x
perc_dxdz = abs((dx-dz)/dz);

disp '--'
disp(['dx = ', num2str(dx)]);
disp(['dz = ', num2str(dz)]);

if (perc_dxdz > 0.01)
    disp 'WARNING: your dx is not your dz. Do you wish to continue?'
    yesorno = input('Do you wish to continue? [yes/no]  ', 's');
    if ( strcmp(yesorno,'n') || strcmp(yesorno,'no') )
        error('You chose to terminate')
        return;
    elseif not( strcmp(yesorno,'yes') || strcmp(yesorno,'y') )
        disp 'you didn''t say yes but we just assume you did anyway.'
        pause(0.3);
    end
end

%% allowed stf frequency

npw = 10; % n samples per wavelength

% minimum wavelength allowed
min_wavelength = min(npw*dz, npw*dz);
max_freq = min(vs(:)) / min_wavelength;

disp(['Your maximum allowed stf frequency is ', num2str(max_freq), ' Hz']);

%% absbound width -- 1x

% typical_wavelength = freq where max of frequency specrum is located

disp(['absbound: ',num2str(width), 'm;  min wavelength: ', num2str(min_wavelength), 'm']);


if (width < min_wavelength) 
    disp('WARNING: absbound less than the minimum allowed wavelength. ');
    disp(['absbound: ',num2str(width), 'm;  min wavelength: ', num2str(min_wavelength),'m']);
    yesorno = input('Do you wish to continue? [yes/no]  ', 's');
    if ( strcmp(yesorno,'n') || strcmp(yesorno,'no') )
        error('You chose to terminate')
%         return;
    elseif not( strcmp(yesorno,'yes') || strcmp(yesorno,'y') )
        disp 'you didn''t say yes but we just assume you did anyway.'
    end
end

if (width > 0.15*min(Lx,Lz))
    disp('WARNING: absbound more than 15% of domain size. ');
    disp(['absbound: ',num2str(width), 'm; Lx: ', num2str(Lx), 'm; Lz: ', num2str(Lz)]);
%     yesorno = input('Do you wish to continue? [yes/no]  ', 's');
%     if ( strcmp(yesorno,'n') || strcmp(yesorno,'no') )
%         error('You chose to terminate')
% %         return;
%     elseif not( strcmp(yesorno,'yes') || strcmp(yesorno,'y') )
%         disp 'you didn''t say yes but we just assume you did anyway.'
%     end
end
    
%% signal reaches receiver?

% may need to change this if I make ONE slow element in the box.
% e.g. calculate the 90 percentile of vp / vs or something.

totaltime =   nt*dt;
dist = sqrt((rec_x-src_x).^2+(rec_z-src_z).^2);

timeneeded_vp =  max(dist(:)) / min(vp(:));
timeneeded_vs =  max(dist(:)) / min(vs(:));

if (timeneeded_vp > totaltime)
    disp 'WARNING: Your P signal may not reach a receiver.'
    yesorno = input('Do you wish to continue? [yes/no]  ', 's');
    if ( strcmp(yesorno,'n') || strcmp(yesorno,'no') )
        error('You chose to terminate')
        return;
    elseif not( strcmp(yesorno,'yes') || strcmp(yesorno,'y') )
        disp 'you didn''t say yes but we just assume you did anyway.'
    end
elseif (timeneeded_vs > totaltime)
    disp 'Your S signal may not reach a receiver.'
    yesorno = input('Do you wish to continue? [yes/no]  ', 's');
    if ( strcmp(yesorno,'n') || strcmp(yesorno,'no') )
        error('You chose to terminate')
        return;
    elseif not( strcmp(yesorno,'yes') || strcmp(yesorno,'y') )
        disp 'you didn''t say yes but we just assume you did anyway.'
    end
end