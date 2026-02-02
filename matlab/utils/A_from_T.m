function A = A_from_T(T)
%A_FROM_T Compute ellipsoid matrix A = T^{-T} T^{-1} given upper triangular T.
X = T \ eye(3);
A = X.'*X;
A = 0.5*(A + A.');
end
