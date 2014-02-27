function [X,Z,K_rho]=run_noise_structure_kernel

%==========================================================================
% run simulation to compute sensitivity kernel for rho (only one-sided)
%
% output:
%--------
% X, Z: coordinate axes
% K_rho: sensitivity kernel
%
%==========================================================================

%==========================================================================
% set paths and read input
%==========================================================================

path(path,'propagation/');
path(path,'../input/');
path(path,'../input/interferometry');

input_parameters;

load cm_velocity;

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

[mu,rho]=define_material_parameters(nx,nz,model_type);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
 
plot_model;

%- time axis --------------------------------------------------------------
    
t=-(nt-1)*dt:dt:(nt-1)*dt;
    
%- read adjoint source locations ------------------------------------------

fid=fopen([adjoint_source_path 'source_locations'],'r');
adsrc_x=zeros(1);
adsrc_z=zeros(1);

k=1;
while (feof(fid)==0)
    adsrc_x(k)=fscanf(fid,'%g',1);
    adsrc_z(k)=fscanf(fid,'%g',1);
    fgetl(fid);
    k=k+1;
end

fclose(fid);

%- read adjoint source time functions and reverse time axis ---------------

nt=length(t);
ns=length(adsrc_x);
stf=zeros(ns,nt);

for n=1:ns
    fid=fopen([adjoint_source_path '/src_' num2str(n)],'r');
    stf(n,nt:-1:1)=fscanf(fid,'%g',nt);
end
    
%- compute indices for adjoint source locations ---------------------------
    
adsrc_x_id=zeros(1,ns);
adsrc_z_id=zeros(1,ns);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:ns

    adsrc_x_id(i)=min(find(min(abs(x-adsrc_x(i)))==abs(x-adsrc_x(i))));
    adsrc_z_id(i)=min(find(min(abs(z-adsrc_z(i)))==abs(z-adsrc_z(i))));

end

%- initialise interferometry ----------------------------------------------
        
input_interferometry;
w_sample=2*pi*f_sample;
dw=w_sample(2)-w_sample(1);

G_1=zeros(nx,nz,length(f_sample));
K_rho=zeros(nx,nz);
             
%- dynamic fields and absorbing boundary field ----------------------------

v=zeros(nx,nz);
absbound=ones(nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

%- initialise absorbing boundary taper a la Cerjan ------------------------

init_absbound;

%==========================================================================
% iterate
%==========================================================================

figure;
set(gca,'FontSize',20);

for n=1:length(t)
    
    %- compute divergence of current stress tensor ------------------------
    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    %- add point sources --------------------------------------------------
    
    for i=1:ns
        DS(adsrc_x_id(i),adsrc_z_id(i))=DS(adsrc_x_id(i),adsrc_z_id(i))+stf(i,n);
    end
    
    %- update velocity field ----------------------------------------------
    
    v=v+dt*DS./rho;
    
    %- apply absorbing boundary taper -------------------------------------
    
    v=v.*absbound;
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    sxy=sxy+dt*mu(1:nx-1,1:nz).*dx_v(v,dx,dz,nx,nz,order);
    szy=szy+dt*mu(:,1:nz-1).*dz_v(v,dx,dz,nx,nz,order);
     
    %- accumulate Fourier transform of the velocity field -----------------
    
    i=sqrt(-1);
    for k=1:length(w_sample)   
        G_1(:,:,k)=G_1(:,:,k)+v(:,:)*exp(-i*w_sample(k)*t(n))*dt;
    end
    
    %- plot velocity field every 4th time step ----------------------------
    
    plot_velocity_field;

end

%==========================================================================
% compute structure kernels 
%==========================================================================

%- load Fourier transformed correlation velocity field
load('../output/interferometry/C_2.mat');

%- accumulate kernel by looping over frequency
for k=1:length(w_sample)
    K_rho=K_rho-G_1(:,:,k).*C_2(:,:,k)*dw;
end

K_rho=real(K_rho);

%==========================================================================
% output 
%==========================================================================

%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    writerObj=VideoWriter(movie_file,'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end