function check_frequency_spectrum(stf,dt)

% check frequency content of stf
%
% input:
% -----
% stf:  the source time function
% dt:   the time step at which it is sampled
%
% output:
% a plot wiht the frequency content of the source time function
%==========================================================================

NFFT = 2^nextpow2(length(stf));
Y = fft(stf,NFFT)/length(stf);
Fs = 1/dt;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(Y(1:NFFT/2+1)));

end