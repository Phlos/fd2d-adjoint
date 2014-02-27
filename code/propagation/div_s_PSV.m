function [outx,outz]=div_s_PSV(sxx,szz,sxz,dx,dz,nx,nz,order)

%==========================================================================
% compute the divergence of the stress tensor for P-SV waves in a medium
% varying only in the X and Z directions (effectively 2D)
% -- Nienke Blom 7 February 2014
%
% input:
%-------
% sxy, szy: stress tensor components
% dx, dz: grid spacings in x- and z-directions
% nx, nz: number of grid points in x- and z-directions
% order: finite-difference order (2 or 4)
%
% output:
%--------
% stress tensor divergence
%==========================================================================

% the strain divergence field for P-SV waves (hence (2,nx,nz) not (3,nx,nz)
% I might include SH but that will only result in more computational effort
% which I really can't be bothered about at the moment. (+ I don't need it)
outx=zeros(nx,nz);
outz=zeros(nx,nz);

if (order==2)

    % DS(1)--------------------------------------------
    % dsxx/dx
    for i=2:nx-1
        outx(i,:)=(sxx(i,:)-sxx(i-1,:))/dx;
    end
    % dsxz/dz
    for j=2:nz-1
        outx(:,j)=outx(:,j)+(sxz(:,j)-sxz(:,j-1))/dz;
    end
    
    % DS(2)--------------------------------------------
    % dsxz/dx
    for i=2:nx-1
        outz(i,:)=(sxz(i,:)-sxz(i-1,:))/dx;
    end
    % dszz/dz
    for j=2:nz-1
        outz(:,j)=outz(:,j)+(szz(:,j)-szz(:,j-1))/dz;
    end
    
    
elseif (order==4)
    
    % DS(1)--------------------------------------------
    % dsxx/dx
    for i=3:nx-2
       outx(i,:)=9*(sxx(i,:)-sxx(i-1,:))/(8*dx) - (sxx(i+1,:)-sxx(i-2,:))/(24*dx);
    end
    % dsxz/dz
    for j=3:nz-2
       outx(:,j)=outx(:,j)+9*(sxz(:,j)-sxz(:,j-1))/(8*dz) - (sxz(:,j+1)-sxz(:,j-2))/(24*dz);
    end
    
    % DS(2)--------------------------------------------
    % dsxz/dx
    for i=3:nx-2
       outz(i,:)=9*(sxz(i,:)-sxz(i-1,:))/(8*dx) - (sxz(i+1,:)-sxz(i-2,:))/(24*dx);
    end
    % dszz/dz
    for j=3:nz-2
       outz(:,j)=outz(:,j)+9*(szz(:,j)-szz(:,j-1))/(8*dz) - (szz(:,j+1)-szz(:,j-2))/(24*dz);
    end
        
end