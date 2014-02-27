function out=div_s(sxy,szy,dx,dz,nx,nz,order)

%==========================================================================
% compute the divergence of the stress tensor
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

out=zeros(nx,nz);

if (order==2)
    % sxy
    for i=2:nx-1
        out(i,:)=(sxy(i,:)-sxy(i-1,:))/dx;
    end
    % szy
    for j=2:nz-1
        out(:,j)=out(:,j)+(szy(:,j)-szy(:,j-1))/dz;
    end
    
elseif (order==4)
    % sxy
    for i=3:nx-2
       out(i,:)=9*(sxy(i,:)-sxy(i-1,:))/(8*dx)-(sxy(i+1,:)-sxy(i-2,:))/(24*dx);
    end
    % szy
    for j=3:nz-2
        out(:,j)=out(:,j)+9*(szy(:,j)-szy(:,j-1))/(8*dz)-(szy(:,j+1)-szy(:,j-2))/(24*dz);
    end
        
end