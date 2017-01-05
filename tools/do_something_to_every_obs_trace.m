function [sObsPerFreq_new] = do_something_to_every_obs_trace(sObsPerFreq) %, SNR)

% prevmax = 0;

butterworth_np = 5;
f_minlist = [0.00667 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667 0.00667];
f_maxlist = [0.00667 0.00833 0.01042 0.01302 0.01628 0.02035 0.02543 0.03179];
sObsPerFreq_new = sObsPerFreq;

% nFreq = numel(sObsPerFreq);
for ii = 1:nFreq
% for ii = 1:numel(f_minlist);
%     freqmax_prev_x = 0; freqmax_prev_z = 0;
    nEvts = numel(sObsPerFreq(1).sEventObs);
    for jj = 1:nEvts
        nRec = numel(sObsPerFreq(1).sEventObs(jj).vel);
        for kk = 1:nRec
            iks = sObsPerFreq(ii).sEventObs(jj).vel{kk}.x;
            zet = sObsPerFreq(ii).sEventObs(jj).vel{kk}.z;
            t_obs = sObsPerFreq(ii).sEventObs(jj).t;
            
            
            %% the "something"  that is done
            
%             % finding the maximum amplitude
%             maxamp = max([abs(iks), abs(zet), prevmax]);
%             if max(abs(iks)) > prevmax
%                 maxtrace = ['Maxtrace: freq ',num2str(ii), ', evt ', num2str(jj), ', rec ', num2str(kk), '.x'];
%             elseif max(abs(zet)) > prevmax
%                 maxtrace = ['Maxtrace: freq ',num2str(ii), ', evt ', num2str(jj), ', rec ', num2str(kk), '.z'];
%             end
%             prevmax = maxamp;
%             % per frequency:
%             freqmax_iks(ii) = max([abs(iks), freqmax_prev_x]);
%             freqmax_zet(ii) = max([abs(zet), freqmax_prev_z]);
%             if max(abs(iks)) > freqmax_prev_x
%                 maxtrace_freqx = ['Maxtrace x:  evt ', num2str(jj), ', rec ', num2str(kk), '.x: ',num2str(freqmax_iks(ii))];
%             elseif max(abs(zet)) > freqmax_prev_z
%                 maxtrace_freqz = ['Maxtrace z:  evt ', num2str(jj), ', rec ', num2str(kk), '.z: ',num2str(freqmax_zet(ii))];
%             end
%             freqmax_prev_x = freqmax_iks(ii);
%             freqmax_prev_z = freqmax_zet(ii);
            
%             % adding noise
%             butterworth_np = 5;
%             
%             % frequencies
%             fmin = sObsPerFreq(ii).f_min;
%             fmax = sObsPerFreq(ii).f_max;
%             
%             % create noise
%             noise_x = randn(size(iks)) - 0.5;
%             noise_z = randn(size(zet)) - 0.5;
%             % filter noise
%             noise_fx = butterworth_lp(noise_x, t_obs, butterworth_np, fmax, 'silent');
%             noise_fx = butterworth_hp(noise_fx, t_obs, butterworth_np, fmin, 'silent');
%             noise_fz = butterworth_lp(noise_z, t_obs, butterworth_np, fmax, 'silent');
%             noise_fz = butterworth_hp(noise_fz, t_obs, butterworth_np, fmin, 'silent');
%             % scale noise
%             noise_fx = noise_fx / max(abs(noise_fx)) * max(abs(iks)) / SNR;
%             noise_fz = noise_fz / max(abs(noise_fz)) * max(abs(zet)) / SNR;
%             % add noise to traces
%             sNew(ii).sEventObs(jj).vel{kk}.x = iks + noise_fx;
%             sNew(ii).sEventObs(jj).vel{kk}.z = zet + noise_fz;

            % filter the traces at the right frequencies
            fmin = f_minlist(ii);
            fmax = f_maxlist(ii);
            iks_f =  butterworth_lp(iks, t_obs, butterworth_np, fmax, 'silent');
            iks_f = butterworth_hp(iks_f, t_obs, butterworth_np, fmin, 'silent');
            zet_f =  butterworth_lp(zet, t_obs, butterworth_np, fmax, 'silent');
            zet_f = butterworth_hp(zet_f, t_obs, butterworth_np, fmin, 'silent');
            
            sObsPerFreq_new(ii).sEventObs(jj).vel{kk}.x = iks_f;
            sObsPerFreq_new(ii).sEventObs(jj).vel{kk}.z = zet_f;

        end
        
    end
%     disp(['Freq ', num2str(ii)]);
%     disp(maxtrace_freqx);
%     disp(maxtrace_freqz);
end

% disp(maxtrace);

end