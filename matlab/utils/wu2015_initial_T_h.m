function [T0, h0, Ae, be, ce] = wu2015_initial_T_h(Y)
%WU2015_INITIAL_T_H Wu & Shi 2015 ellipsoid-based initialization.
%
% Fits an ellipsoid in algebraic form:
%   y^T Ae y + be^T y + ce = 0
% and rescales it so that
%   (y-h0)^T A (y-h0) = 1
% with A SPD. Then A = R^T R and T0 = R^{-1}.
%
% Inputs:
%   Y 3xN samples
%
% Outputs:
%   T0 upper triangular
%   h0 bias
%   Ae,be,ce unscaled algebraic coefficients

x = Y(1,:).';
y = Y(2,:).';
z = Y(3,:).';
N = size(Y,2);

D = [ x.^2, y.^2, z.^2, 2*x.*y, 2*x.*z, 2*y.*z, x, y, z, ones(N,1) ];
[ V, Eval ] = eig(D.'*D);
[~, idx] = min(diag(Eval));
pv = V(:,idx);

a11 = pv(1); a22 = pv(2); a33 = pv(3);
a12 = pv(4); a13 = pv(5); a23 = pv(6);
be = pv(7:9);
ce = pv(10);

Ae = [a11, a12, a13;
      a12, a22, a23;
      a13, a23, a33];
Ae = 0.5*(Ae+Ae.');

% Wu scaling
Den = (be.' * (Ae \ be)) - 4*ce;
alpha = 4 / Den;

A = alpha * Ae;
b = alpha * be;

% Make SPD (sign ambiguity in the implicit quadratic)
if min(eig(A)) <= 0
    A = -A;
    b = -b;
end

h0 = -0.5*(A \ b);
A = 0.5*(A + A.');
R0 = chol(A);           % upper triangular, R0'*R0 = A
T0 = inv(R0);
end
