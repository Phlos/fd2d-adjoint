function Kf = filter_kernels(K,sigma)

% filters a 2D kernel K with a gaussian filter of width sigma. This is more
% efficient than using a regular 2D filter, because instead of going over
% nx*ny points with the filter, one goes once in the x direction, and once
% in the y direction, thus nx+ny operations.
%
% INPUT:
% - K:      Kernel to be filtered (2D array)
% - sigma:  width of Gaussian filter - unit of pixels
%
% OUTPUT:
% - Kf:     Gaussian filtered kernel.
%
% [ sent by Andreas Fichtner on 31-10-2014 ]
% [ adapted on 4-11-2014 (Nienke Blom): 
%   - added intermediate kernel Kfy 
%   - removed grid dependencies            ]


%- PREPARATION

nx=size(K,2);
ny=size(K,1);

yy = [0:ny-1];
xx = [0:nx-1];

Kfy = zeros(size(K));
Kf = zeros(size(K));


%- ACTUAL FILTERING

%- Filtering in y-direction.

for j=1:ny
    
    %- Gaussian.
    g=(1.0/(sqrt(2.0*pi)*sigma)) * exp(-(yy-yy(j)).^2 / (2.0*sigma^2));

    %- Area beneath the Gaussian for normalisation.
    N=sum(g);
    
    %-Compute the convolution.
    for i=1:nx 
        Kfy(j,i)=sum(K(:,i).*g')/N;
    end

end

%- Filtering in x-direction.

for i=1:nx
    
    %- Gaussian.
    g=(1.0/(sqrt(2.0*pi)*sigma)) * exp(-(xx-xx(i)).^2 / (2.0*sigma^2));

    %- Area beneath the Gaussian for normalisation.
    N=sum(g);
    
    %-Compute the convolution.
    for j=1:ny 
        Kf(j,i)=sum(Kfy(j,:).*g)/N;
    end

end
    
end

