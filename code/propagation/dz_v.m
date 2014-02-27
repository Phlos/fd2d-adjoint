function out=dz_v(v,dx,dz,nx,nz,order)

%==========================================================================
% compute z-derivative of velocity
%
% input:
%-------
% v: velocity
% dx, dz: grid spacings in x- and z-directions
% nx, nz: number of grid points in x- and z-directions
% order: finite-difference order (2, 4, 6 or 8)
%
% output:
%--------
% velocity derivative in z-directions
%==========================================================================

out=zeros(nx,nz-1);

if (order==2)

    for j=1:nz-1
        out(:,j)=(v(:,j+1)-v(:,j))/dz;
    end
    
elseif (order==4)
    
    for j=2:nz-2
        out(:,j)=9*(v(:,j+1)-v(:,j))/(8*dz)-(v(:,j+2)-v(:,j-1))/(24*dz);
    end
    
end