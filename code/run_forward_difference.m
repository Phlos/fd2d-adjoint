function run_forward_difference(plot_mode,v_max)

%==========================================================================
% run forward simulation for regular and perturbed wave field 
%
% input:
%-------
% plot_mode: show wavefield simulation when plot_mode=='plot'
% v_max: maximum values of the colorscale for plotting
%
%==========================================================================

%==========================================================================
% read input and set paths
%==========================================================================

input_parameters;

path(path,'helper_programmes/');

if strcmp(plot_mode,'plot')
    h_vel=figure;
    load cm;
end

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

%- unperturbed medium
[mu,rho]=define_material_parameters(nx,nz,1);

%- perturbed medium
[mu_pert,rho_pert]=define_material_parameters(nx,nz,2);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

if strcmp(plot_mode,'plot')
    
    h_model=figure;
    
    subplot(2,2,1)
    pcolor(X,Z,mu');
    axis image
    shading flat
    title('mu [N/m^2]')
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    
    subplot(2,2,2)
    pcolor(X,Z,rho');
    axis image
    shading flat
    title('rho [kg/m^3]')
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    
    subplot(2,2,3)
    pcolor(X,Z,mu_pert');
    axis image
    shading flat
    title('mu perturbed [N/m^2]')
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
    
    subplot(2,2,4)
    pcolor(X,Z,rho_pert');
    axis image
    shading flat
    title('rho perturbed [kg/m^3]')
    xlabel('x [m]');
    ylabel('z [m]');
    colorbar
        
end
    
%- dynamic fields ---------------------------------------------------------

v=zeros(nx,nz);
v_pert=zeros(nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

sxy_pert=zeros(nx-1,nz);
szy_pert=zeros(nx,nz-1);

%- read seismic source locations ------------------------------------------

fid=fopen([source_path 'source_locations'],'r');
src_x=zeros(1);
src_z=zeros(1);

k=1;
while (feof(fid)==0)
    src_x(k)=fscanf(fid,'%g',1);
    src_z(k)=fscanf(fid,'%g',1);
    fgetl(fid);
    k=k+1;
end

fclose(fid);

%- read source time functions ---------------------------------------------

ns=length(src_x);
stf=zeros(ns,nt);

for n=1:ns
    fid=fopen([source_path 'src_' num2str(n)],'r');
    stf(n,1:nt)=fscanf(fid,'%g',nt);
end

%- compute indices for source locations -----------------------------------

src_x_id=zeros(1,ns);
src_z_id=zeros(1,ns);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:ns

    src_x_id(i)=min(find(min(abs(x-src_x(i)))==abs(x-src_x(i))));
    src_z_id(i)=min(find(min(abs(z-src_z(i)))==abs(z-src_z(i))));

end

%==========================================================================
% iterate
%==========================================================================

t=0;

figure(h_vel);

for n=1:nt
    
    %- compute divergence of current stress tensor and add external forces
    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    DS_pert=div_s(sxy_pert,szy_pert,dx,dz,nx,nz,order);
    
    for i=1:ns
        DS(src_x_id(i),src_z_id(i))=DS(src_x_id(i),src_z_id(i))+stf(i,n);
        DS_pert(src_x_id(i),src_z_id(i))=DS_pert(src_x_id(i),src_z_id(i))+stf(i,n);
    end
    
    %- update velocity field ----------------------------------------------
    
    v=v+dt*DS./rho;
    v_pert=v_pert+dt*DS_pert./rho_pert;
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    sxy=sxy+dt*mu(1:nx-1,1:nz).*dx_v(v,dx,dz,nx,nz,order);
    szy=szy+dt*mu(:,1:nz-1).*dz_v(v,dx,dz,nx,nz,order);
    
    sxy_pert=sxy_pert+dt*mu_pert(1:nx-1,1:nz).*dx_v(v_pert,dx,dz,nx,nz,order);
    szy_pert=szy_pert+dt*mu_pert(:,1:nz-1).*dz_v(v_pert,dx,dz,nx,nz,order);
    
    %- update time --------------------------------------------------------
    
    t=t+dt;
    
    %- plot velocity field ------------------------------------------------
    
    if strcmp(plot_mode,'plot')
    
        subplot(1,3,1)
        pcolor(X,Z,v');
    
        caxis([-v_max,v_max]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('velocity field [m/s]');
        
        subplot(1,3,2)
        pcolor(X,Z,v_pert');
    
        caxis([-v_max,v_max]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('perturbed velocity field [m/s]');
        
        subplot(1,3,3)
        pcolor(X,Z,v_pert'-v');
    
        caxis([-v_max/2,v_max/2]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('perturbed - original velocity field [m/s]');
        
        %colorbar
        
        hold on
        for k=1:ns
            plot(src_x(k),src_z(k),'kx')
        end
        hold off

        pause(0.02)
        
    end
    
end