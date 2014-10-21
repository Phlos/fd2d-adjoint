function C = add_crop_matrix(A, B, dx, dz)

% C = add_crop_matrix(A,B, dx,dz)
% function that adds matrix B to matrix A, with the center of B at 
% (dx, dz) inside A. B is cropped where it falls outside of A.

% A = zeros(5,8)
% B = ones(3,3)

[b.width, b.height] = size(B);
[a.width, a.height] = size(A);

if rem(b.width,2) == 0
    dx = dx + 0.5;
end
if rem(b.height,2) == 0
    dz = dz + 0.5;
end
    
% axstartnoround = dx - 0.5* b.width
% x start & stop values
a.xstart = dx - 0.5* b.width + 0.5; % + 0.5 
% axstart = a.xstart
% halfbwidth = 0.5*b.width
% rounddxminushalfbwid = round(dx - 0.5* b.width)
if a.xstart < 1
    b.xstart = abs(a.xstart) + 2;
    a.xstart = 1;
else
    b.xstart = 1;
end
% a.xstart
% b.xstart
a.xstop = dx + 0.5*b.width - 0.5;
% axstop = a.xstop
if a.xstop > a.width
    b.xstop = b.width - (a.xstop - a.width);
    a.xstop = a.width;
else
    b.xstop = b.width;
end
% a.xstop
% b.xstop
% z start & stop values
a.zstart = dz - 0.5*b.height + 0.5;
if a.zstart < 1
    b.zstart = abs(a.zstart) + 2;
    a.zstart = 1;
else
    b.zstart = 1;
end
a.zstop = dz + 0.5*b.height - 0.5;
if a.zstop > a.height
    b.zstop = b.height - (a.zstop - a.height);
    a.zstop = a.height;
else
    b.zstop = b.height;
end

% a
% b

C = A;
C(a.xstart:a.xstop , a.zstart:a.zstop) = C(a.xstart:a.xstop , a.zstart:a.zstop) ...
                                       + B(b.xstart:b.xstop , b.zstart:b.zstop);

end