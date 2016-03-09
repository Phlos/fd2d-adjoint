function [sObs_noised] = add_noise_to_traces(sObsPerFreq, SNR)

    % function to add noise to the (observed) data.
    % [sObs_noised] = add_noise_to_traces(sObsPerFreq, SNR);
    % INPUT:
    % sObsPerFreq: struct with sObsPerFrec(ii).sEventObs(jj).vel{kk}.x
    % SNR:         signal-to-noise ratio (max amplitude of signal / max
    % amplitude of noise)
    %
    % OUTPUT:
    % sObs_noised: same struct as sObsPerFreq, but with noise
    %
    % -- N.A.Blom, 24-2-2016
    
    %% preparation
    
    % filtering
    butterworth_np = 5;
    
    sObs_noised = sObsPerFreq;
    
    %% create noise traces per source-receiver pair
    % (taken from frequency 1)
    
    disp('Creating noise traces...');
    nEvts = numel(sObsPerFreq(1).sEventObs);
    for jj = 1:nEvts
        nRec = numel(sObsPerFreq(1).sEventObs(jj).vel);
        for kk = 1:nRec
            
            % traces
            iks = sObsPerFreq(1).sEventObs(jj).vel{kk}.x;
            zet = sObsPerFreq(1).sEventObs(jj).vel{kk}.z;
            % create noise traces per frequency
            noise.evt(jj).rec(kk).x = randn(size(iks)) - 0.5;
            noise.evt(jj).rec(kk).z = randn(size(zet)) - 0.5;
            
        end
    end
    
    %% add (filtered) noise to the traces
    
    disp('Adding noise to OBS data...');
    nFreq = numel(sObsPerFreq);
    for ii = 1:nFreq
        nEvts = numel(sObsPerFreq(ii).sEventObs);
        for jj = 1:nEvts
            nRec = numel(sObsPerFreq(ii).sEventObs(jj).vel);
            for kk = 1:nRec
                
                iks = sObsPerFreq(ii).sEventObs(jj).vel{kk}.x;
                zet = sObsPerFreq(ii).sEventObs(jj).vel{kk}.z;
                t_obs = sObsPerFreq(ii).sEventObs(jj).t;
                
                
                %% adding noise
                
                % frequencies
                fmin = sObsPerFreq(ii).f_min;
                fmax = sObsPerFreq(ii).f_max;
                
                % filter noise
                noise_x = noise.evt(jj).rec(kk).x;
                noise_fx = butterworth_lp(noise_x, t_obs, butterworth_np, fmax, 'silent');
                noise_fx = butterworth_hp(noise_fx, t_obs, butterworth_np, fmin, 'silent');
                noise_z = noise.evt(jj).rec(kk).z;
                noise_fz = butterworth_lp(noise_z, t_obs, butterworth_np, fmax, 'silent');
                noise_fz = butterworth_hp(noise_fz, t_obs, butterworth_np, fmin, 'silent');
                % scale noise
                noise_fx = noise_fx / max(abs(noise_fx)) * max(abs(iks)) / SNR;
                noise_fz = noise_fz / max(abs(noise_fz)) * max(abs(zet)) / SNR;
                % add noise to traces
                sObs_noised(ii).sEventObs(jj).vel{kk}.x = iks + noise_fx;
                sNew(ii).sEventObs(jj).vel{kk}.z = zet + noise_fz;
                
            end
        end
    end
%     
%     for ii = 1:nFreq; 
%         sObs_noised(ii).sEventObs = sObs_noised(ii).sEventObs; end;
    
end