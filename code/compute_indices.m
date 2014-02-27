% calculates the indices at which objects (sources, receivers) reside,
% based on their location in true coordinates

function [coord_x_id, coord_z_id, nthings] = compute_indices(coord_x,coord_z,Lx,Lz,dx,dz)

    nthings=length(coord_x);
    
    coord_x_id=zeros(1,nthings);
    coord_z_id=zeros(1,nthings);

    x=0:dx:Lx;
    z=0:dz:Lz;

    for i=1:nthings

        coord_x_id(i)=min(find(min(abs(x-coord_x(i)))==abs(x-coord_x(i))));
        coord_z_id(i)=min(find(min(abs(z-coord_z(i)))==abs(z-coord_z(i))));

    end

end