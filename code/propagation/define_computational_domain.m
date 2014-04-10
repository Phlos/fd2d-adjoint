function [X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz)

%==========================================================================
% define the geometry and discretisation of the computational domain
%
% input:
%-------
% Lx, Lz: extensions of the computational doman [m] in x- and z-directions
% nx, nz: number of grid points in x- and z-directions
%
% output:
%--------
% X, Z: coordinae matrices
% dx, dz: grid sizes in x- and z-directions
%==========================================================================

dx=Lx/(nx-1);
dz=Lz/(nz-1);

[X,Z]=meshgrid(0:dx:Lx,0:dz:Lz);

