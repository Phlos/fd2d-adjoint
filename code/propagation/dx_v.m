function out=dx_v(v,dx,dz,nx,nz,order)

%==========================================================================
% compute x-derivative of velocity
%
% input:
%-------
% v: velocity
% dx, dz: grid spacings in x- and z-directions
% nx, nz: number of grid points in x- and z-directions
% order: finite-difference order (2, 4, 6 or 8) ---> welll you only
%        implemented the orders up to order 4 my friend...
%
% output:
%--------
% velocity derivative in x-directions
%==========================================================================

out=zeros(nx,nz);

if (order==2)

    for i=1:nx-1
        out(i,:)=(v(i+1,:)-v(i,:))/dx;
    end
    
elseif (order==4)
   
    for i=2:nx-2
        out(i,:)=9*(v(i+1,:)-v(i,:))/(8*dx)-(v(i+2,:)-v(i-1,:))/(24*dx);
    end
    
end