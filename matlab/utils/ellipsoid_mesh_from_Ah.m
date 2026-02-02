function [Xr, Yr, Zr] = ellipsoid_mesh_from_Ah(A, h, n)
%ELLIPSOID_MESH_FROM_AH Create a rotated ellipsoid mesh from quadratic form.
%
% (y-h)'*A*(y-h) = 1 with A SPD.
%
% Uses MATLAB ellipsoid() to generate an axis-aligned ellipsoid in principal
% coordinates and then rotates it by the eigenvectors of A.
%
% Inputs:
%   A 3x3 SPD matrix
%   h 3x1 center
%   n mesh density for ellipsoid() (default 40)
%
% Outputs:
%   Xr, Yr, Zr arrays suitable for surf()

if nargin < 3 || isempty(n)
    n = 40;
end

A = 0.5*(A + A.');
[V,D] = eig(A);
d = diag(D);
if any(d <= 0)
    error('A must be SPD. Smallest eigenvalue = %g', min(d));
end
r = 1 ./ sqrt(d);

[X,Y,Z] = ellipsoid(0,0,0, r(1), r(2), r(3), n);

pts = V * [X(:)'; Y(:)'; Z(:)'];
Xr = reshape(pts(1,:), size(X)) + h(1);
Yr = reshape(pts(2,:), size(Y)) + h(2);
Zr = reshape(pts(3,:), size(Z)) + h(3);
end
