%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- number of noise sources
n_noise_sources=3;

%- characteristics of the noise spectrum ----------------------------------
%- only needed in this routine --------------------------------------------
f_peak=1.0/16.0;       % peak frequency in Hz
bandwidth=0.35/16.0;    % bandwidth in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_spectrum=zeros(length(f_sample),n_noise_sources);

%- reference noise spectrum for whitening

noise_spectrum_ref=0.45*exp(-(abs(f_sample)-0.07).^2/bandwidth^2)+1.0*exp(-(abs(f_sample)-0.12).^2/(0.5*bandwidth)^2)+0.036*[f_sample<0.12].*[f_sample>0.02];
taper=0.02*[f_sample<0.15].*[f_sample>0.02];

%- spectrum for source region 1 -------------------------------------------

noise_spectrum(:,1)=taper.*(exp(-(abs(f_sample)-0.07).^2/bandwidth^2))./noise_spectrum_ref;

figure
set(gca,'FontSize',20)
hold on

plot(f_sample,noise_spectrum(:,1),'k');
xlabel('frequency [Hz]','FontSize',20);
title('noise power-spectral density for source region 1' ,'FontSize',20);

%- spectrum for source region 2 -------------------------------------------

noise_spectrum(:,2)=taper.*(exp(-(abs(f_sample)-0.12).^2/(0.5*bandwidth)^2))./noise_spectrum_ref;

figure
set(gca,'FontSize',20)
hold on

plot(f_sample,noise_spectrum(:,2),'k');
xlabel('frequency [Hz]','FontSize',20);
title('noise power-spectral density for source region 2' ,'FontSize',20);

%- spectrum for source region 3 -------------------------------------------

noise_spectrum(:,3)=taper.*(0.02*[f_sample<0.12].*[f_sample>0.02])./noise_spectrum_ref;

figure
set(gca,'FontSize',20)
hold on

plot(f_sample,noise_spectrum(:,3),'k');
xlabel('frequency [Hz]','FontSize',20);
title('noise power-spectral density for source region 3' ,'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geographic power-spectral density distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_source_distribution=zeros(nx,nz,n_noise_sources);

%- noise source geography for region 1 ------------------------------------

noise_source_distribution(:,:,1)=(exp(-((X-0.6e6).^2+(Z-0.6e6).^2)/(0.2e6)^2))';

figure;
set(gca,'FontSize',20);
load cm_psd

pcolor(X,Z,noise_source_distribution(:,:,1)');
shading interp
colormap(cm_psd)
xlabel('x [m]','FontSize',20);
ylabel('z [m]','FontSize',20);
title('power-spectral density distribution of noise sources for region 1','FontSize',20);

%- noise source geography for region 2 ------------------------------------

noise_source_distribution(:,:,2)=(exp(-((X-0.4e6).^2+(Z-1.1e6).^2)/(0.3e6)^2))';

figure;
set(gca,'FontSize',20);
load cm_psd

pcolor(X,Z,noise_source_distribution(:,:,2)');
shading interp
colormap(cm_psd)
xlabel('x [m]','FontSize',20);
ylabel('z [m]','FontSize',20);
title('power-spectral density distribution of noise sources for region 2','FontSize',20);


%- noise source geography for region 3 ------------------------------------

noise_source_distribution(:,:,3)=1.0;

figure;
set(gca,'FontSize',20);
load cm_psd

pcolor(X,Z,noise_source_distribution(:,:,3)');
shading interp
colormap(cm_psd)
xlabel('x [m]','FontSize',20);
ylabel('z [m]','FontSize',20);
title('power-spectral density distribution of noise sources for region 3','FontSize',20);

