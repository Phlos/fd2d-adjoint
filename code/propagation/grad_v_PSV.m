function [dvxdx,dvxdz,dvzdx,dvzdz]=grad_v_PSV(vx,vz,dx,dz,nx,nz,order)

%==========================================================================
% compute the gradient of the velocity field for P-SV waves 
% in a medium varying only in the X and Z directions (effectively 2D)
% -- Nienke Blom 7 February 2014
% NOTE: can be used for any 2D vector field, really. Done so for kernel
% calculations in adjoint mode. -- NAB 13 March 2014
%
% input:
%-------
% vx, vz: velocity components in x and z directions
% dx, dz: grid spacings in x- and z-directions
% nx, nz: number of grid points in x- and z-directions
% order: finite-difference order (2 or 4)
%
% output:
%--------
% gradient of velocity field
%==========================================================================

% initialise the velocity gradient tensor. At each grid point 2x2 variables
% are defined: dvx/dx, dvx/dz, dvz/dx, dvz/dz
%  out=zeros(2,2,nx,nz); % GRRRRRR DOES NOT WORK ---> out(i,j,k,l) = dvi/dxj(k,l)
dvxdx=zeros(nx,nz);
dvxdz=zeros(nx,nz);
dvzdx=zeros(nx,nz);
dvzdz=zeros(nx,nz);

if (order==2)  % deze klopt volgens mijnnnnn
    % dvx/dx
    for i=1:nx-1
        dvxdx(i,:)=(vx(i+1,:)-vx(i,:))/dx;
    end
    % dvx/dz
    for j=1:nz-1
        dvxdz(:,j)=(vx(:,j+1)-vx(:,j))/dz;
    end
    % dvz/dx
    for i=1:nx-1
        dvzdx(i,:)=(vz(i+1,:)-vz(i,:))/dx;
    end
    % dvz/dz
    for j=1:nz-1
        dvzdz(:,j)=(vz(:,j+1)-vz(:,j))/dz;
    end
    
    dx_v
elseif (order==4)   % Nou deze klopt ook... denk ik... ik kan geen i en j meer zien.aarhghghghghghhgghhghghgfhg
    % dvx/dx
    for i=2:nx-2
        dvxdx(i,:)=9*(vx(i+1,:)-vx(i,:))/(8*dx)-(vx(i+2,:)-vx(i-1,:))/(24*dx);
    end
    % dvx/dz
    for j=2:nz-2
        dvxdz(:,j)=9*(vx(:,j+1)-vx(:,j))/(8*dz)-(vx(:,j+2)-vx(:,j-1))/(24*dz);
    end
    % dvz/dx
    for i=2:nx-2
        dvzdx(i,:)=9*(vz(i+1,:)-vz(i,:))/(8*dx)-(vz(i+2,:)-vz(i-1,:))/(24*dx);
    end
    % dvx/dz
    for j=2:nz-2
        dvzdz(:,j)=9*(vz(:,j+1)-vz(:,j))/(8*dz)-(vz(:,j+2)-vz(:,j-1))/(24*dz);
    end
    
end
    
    
