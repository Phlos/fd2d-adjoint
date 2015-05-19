function [adstf, misfit] = calc_misfitseis_adstf(misfit_tiep, t, v_rec, varargin)

% function calculating the seismic misfit and the adjoint stf
%
% SYNTAX:
% [adstf, misfit] = calc_misfitseis_adstf(misfit_type, t, v_rec)
%                   automatic picking of whole trace (minus taper) 
%                   v_obs is assumed zero, picking_mode is 'auto'
% [adstf, misfit] = calc_misfitseis_adstf(misfit_type, t, v_rec, v_obs)
%                   automatic picking of whole trace (minus taper)
%                   v_obs is supplied (necessary for waveform_difference),
%                   picking_mode is assumed 'auto'
% [adstf, misfit] = calc_misfitseis_adstf(misfit_type, t, picking_mode)
%                   v_obs is assumed zero, picking_mode supplied
% [adstf, misfit] = calc_misfitseis_adstf(misfit_type, t, v_obs, picking_mode)
%                   v_obs is supplied, picking_mode supplied
% 
% INPUT:
% - misfit_type:    which misfit functional you're using
%                   waveform_difference:    L2 norm of velocity seismograms
%                   cc_time_shift:          cross-correlation time shift
%                   cc_time_shift_x:        x component picking
%                   cc_time_shift_z:        z component picking
% - t:              time axis, needed if displacement needs to be calc'd 
%                   from vel or the other way around
% - v_obs:          observed seismograms (not always needed)
%                   cell array of shape v_obs{irec}.componentxyz(t);
% - v_rec:          seismograms obtained from wave propagation in current 
%                   model
%                   cell array of shape v_rec{irec}.componentxyz(t);
% - picking_mode:   'auto' or 'manual' (if omitted, will be 'auto')
%
% OUTPUT:
% - adstf:          adstf{irec}.componentxyz(t).
% - misfit:         struct misfit.total misfit.componentxyz
%
% -- Nienke Blom, 28 January 2015


%- prepare 
input_parameters; % should get nrec here, also nt and dt
[~, ~, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);

switch misfit_tiep
    case 'cc_time_shift_x'
        misfit_type = 'cc_time_shift';
        compwel = 'x';
    case 'cc_time_shift_z'
        misfit_type = 'cc_time_shift';
        compwel = 'z';
    otherwise
        misfit_type = misfit_tiep;
        compwel = {'x', 'y', 'z'};
end

switch misfit_type
    case 'waveform_difference'
%         disp 'we''re doing waveform differences'
    case 'cc_time_shift'
%         disp 'we''re doing cc time shifts'
    otherwise
        error('unknown misfit type!');
end

%- get v_obs from varargin
[vobspresent, v_obs, picking_mode] = checkargs(varargin(:));

%- prepare figure for manual if manual picking is selected
if strcmp(picking_mode, 'manual')
    fig_manual = figure;
else
    fig_manual = NaN;
end



%- loop over receivers
nrec = size(v_rec,2);
for irec = 1:nrec
    comp = fieldnames(v_rec{irec});
    for icomp = 1:length(comp)
        if find(strcmp(comp{icomp},compwel))
            % make temporary 1d arrays of currently used recording
            vrec = v_rec{irec}.(comp{icomp});
            if strcmp(vobspresent,'yes')
                vobs = v_obs{irec}.(comp{icomp});
            else
                vobs = zeros(size(t));
            end
            
            %- pick left and right manual or auto
            if strcmp(picking_mode, 'manual')
                [left, right] = pick_adstf_manual(t, vrec, vobs, irec, fig_manual,misfit_type);
                disp(['left: ',num2str(left),'   right: ',num2str(right)]);
            else
                left = t(1);
                right = t(end);
            end
            
            % %- taper vrec (and vobs)
            % no, v_rec & v_obs should not be tapered, since that'll result
            % in weirdness when, after tapering, the recordings are changed
            % from velocity to displacement (needed for example @ waveform
            % difference misfit functionals)
            
    %         bips = figure;
    %         subplot(3,1,1); plot(t,vrec,t,vobs);
            
            %- the actual adstf creation
            
            % first get the adstf_nonreversed
            if strcmp(misfit_type,'waveform_difference')
                [misfit_n,adstf_temp] = misfit_wavef_L2(vrec,vobs,t);
            elseif strcmp(misfit_type,'cc_time_shift')
                [misfit_n,adstf_temp] = misfit_cc_tshft(vrec,vobs,t);
            end
            
    %         figure(bips)
    %         subplot(3,1,2); plot(t,adstf_temp); title('adstf');
            
            % taper adstf_temp
            taper_width=t(end)/10;   % why divide by ten??
            lmax = taper_width;
            rmax = t(end) - taper_width;
            if (left < taper_width) || (right > t(end) - taper_width)
                left = max(left,lmax);
                right = min(right,rmax);
            end

            adstf_temp = taper(fliplr(adstf_temp),t,left,right,taper_width);
            adstf_temp = fliplr(adstf_temp);
    %         figure(bips)
    %         subplot(3,1,3); plot(t,adstf_temp); title('tapered'),
            
            % if manual: plot temp adstf before & after time reversal
            if strcmp(picking_mode, 'manual')
                figure(fig_manual);
                subplot(3,1,2);
                plot(t,fliplr(adstf_temp),'k')
                xlabel('t [s]')
                title(['adjoint source before time reversal'])
                
                subplot(3,1,3);
                plot(t,adstf_temp,'k')
                xlabel('t [s]')
                title(['adjoint source after time reversal'])
                pause(10.5);
            end
            
            % give the adstf the right magnitude so that the spatial integral
            % is 1 (because the stf/adstf are spatial delta functions)
            adstf_temp = adstf_temp / dx / dz;
            
            adstf{irec}.(comp{icomp}) = adstf_temp;
            
            % misfit
            misfit.(comp{icomp})(irec) = misfit_n;
            
        else
            adstf{irec}.(comp{icomp}) = zeros(size(v_rec{irec}.(comp{icomp})));
            misfit.(comp{icomp})(irec) = 0;
            
        end
    end
    

end

    % calculate the total misfit
    comp = fieldnames(misfit);
    misfit.total = 0;
    for icomp = 1:length(comp)
        misfit.total = misfit.total + sum(misfit.(comp{icomp}));
    end



end



function [vobspresent, v_obs, picking_mode] = checkargs(arg)

narg = length(arg);

if narg == 0
    vobspresent = 'no';
    v_obs = NaN;
    picking_mode = 'auto';
    
elseif narg == 1
    if iscell(arg{1})
        vobspresent = 'yes';
        v_obs = arg{1};
        picking_mode = 'auto';
    elseif ischar(arg{1})
        vobspresent = 'no';
        v_obs = NaN;
        picking_mode = arg{1};
    else
        error('the input to calc_misfitseis_adstf was not understood')
    end
    
elseif narg == 2
    vobspresent = 'yes';
    v_obs = arg{1};
    picking_mode = arg{2};
    
else 
    error('the input to calc_misfitseis_adstf was not understood')
end

if ~( strcmp(picking_mode,'manual') || strcmp(picking_mode,'auto') )
    warning(['warning, given picking mode was "',picking_mode, '" -- we assume auto'])
end

end



function [left, right] = pick_adstf_manual(t, v_rec, v_obs, n, fig_manual,misfit_type)

%- preparation
dt = t(2)-t(1);

if strcmp(misfit_type,'waveform_difference')
    u=cumsum(v_rec,2)*dt;
    v_rec=u;
    u_0 = cumsum(v_obs,2)*dt;
    v_obs = u_0;
end

%- plotting and picking

fprintf(1,'station number %d\n',n)

figure(fig_manual);
clf;
subplot(3,1,1)
plot(t,v_rec(:),'k')
hold on
plot(t,v_obs(:),'r')
plot(t,v_rec(:)-v_obs(:),'k--','LineWidth',2)
% hold off

title(['receiver ' num2str(n) ' ,synth - black, obs - red, diff - dashed'])
xlabel('t [s]')
if strcmp(misfit_type,'waveform_difference')
    ylabel('displacement [m]')
else
    ylabel('velocity [m/s]')
end
    

disp('select left window');
[left,~]=ginput(1);
disp('select right window');
[right,~]=ginput(1);

% hold off;
ei_as = get(gca,'ylim');
plot([left, left], ei_as, '--');
plot([right, right], ei_as, '--');
hold off

end



function [misfit,adstf]=misfit_cc_tshft(v,v_0,t)

%- compute cross-correlation time shift -----------------------------------
%  uses velocity seismograms
%
% function [misfit,adsrc]=cc_time_shift(u,u_0,t)
%
% v:    synthetic velocity seismogram
% v_0:  observed velocity seismogram
% t:    time axis

%- initialisations --------------------------------------------------------
dt=t(2)-t(1);
% path(path,'tf/');

%- compute time shift -----------------------------------------------------

if sum(v_0==0)==length(t) % if all v_0 are zero
    T=1.0;
else
%     [cc,t_cc]=cross_correlation(v_0,v,t);
    [cc,t_cc] = xcorr(v_0,v);
    t_cc = t_cc*dt;
%     figure; plot(t_cc,cc);
    T=t_cc(find(max(cc)==cc));
end

%- compute adjoint source time function -----------------------------------


% nt=length(t);
% v=zeros(1,nt);

% v(1:nt-1)=diff(u)/dt;

adstf=T*v/(sum(v.*v)*dt);

%- correction of adjoint source time function for velocity output -
%  (because the signal is time-reversed, +velocity in the positive
%  time direction becomes -velocity in the negative time direction)
adstf = -1 * adstf;

%- time-reverse the adjoint stf
adstf = fliplr(adstf);

%- compute misfit ---------------------------------------------------------

misfit=T^2/2.0;
    
end

%- compute cross correlation function -------------------------------------
%
% function [cc,t_cc]=cross_correlation(f,g,t)
%
% cc_i = sum_j f^*(j) g(j+i)

function [cc,t_cc]=cross_correlation(f,g,t)

%- initialisations --------------------------------------------------------

n=length(f);
cc=zeros(2*n+1);

dt=t(2)-t(1);
t_cc=-n*dt:dt:n*dt;

%- compute brute-force correlation function -------------------------------

for i=-n:n
    for j=1:length(f)
        cc(n+i+1)=cc(n+i+1)+conj(f(j))*g(j+i);
    end
end

end



function [misfit,adstf]=misfit_wavef_L2(v,v_0,t)

%- compute L2 waveform difference -----------------------------------------
% uses displacement seismograms
%
% function [misfit,adsrc]=waveform_difference(u,u_0,t)
%
% v: synthetic velolcity seismogram
% v_0: observed velocity seismogram
% t: time axis

dt = t(2) - t(1);

u = cumsum(v,2)*dt;
u_0 = cumsum(v_0,2)*dt;

adstf=u-u_0;

%- time-reverse the adjoint stf
adstf = fliplr(adstf);

misfit=sum(adstf.*adstf) * dt;


end
