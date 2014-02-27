function plot_noise_source_kernels(X,Z,K_s)

%==========================================================================
%- plot noise source kernels as a function of frequency 
%==========================================================================

%- initialisations --------------------------------------------------------

path(path,'../input/');
load cm_velocity;

input_parameters;
input_interferometry;

%- decide on a frequency --------------------------------------------------

fprintf(1,'f_min=%g Hz, f_max=%g Hz, df=%g Hz\n',min(f_sample),max(f_sample),f_sample(2)-f_sample(1));
f_min=input('min frequency in Hz: ');
f_max=input('max frequency in Hz: ');

ind_f_min=find(min(abs(f_min-f_sample))==abs(f_min-f_sample));
ind_f_max=find(min(abs(f_max-f_sample))==abs(f_max-f_sample));

%- plot source and receiver positions -------------------------------------

figure
hold on
    
for k=1:length(src_x)
    plot(src_x(k),src_z(k),'kx')
end

for k=1:length(rec_x)
    plot(rec_x(k),rec_z(k),'ko')
end
       
%- plot noise source kernel -----------------------------------------------

K_sum=sum(K_s(:,:,ind_f_min:ind_f_max),3);

pcolor(X,Z,K_sum');
set(gca,'FontSize',20);
axis equal

%- axis, scaling, etc. ----------------------------------------------------

colormap(cm);
shading interp
colorbar

m=max(max(K_sum));
caxis([-m m]);

xlabel('x [m]','FontSize',20);
ylabel('z [m]','FontSize',20);


