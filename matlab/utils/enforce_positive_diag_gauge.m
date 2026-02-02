function [Tpos, Mpos, Q] = enforce_positive_diag_gauge(T, M)
%ENFORCE_POSITIVE_DIAG_GAUGE Enforce positive diagonal on an upper triangular T.
% Applies a diagonal reflection Q = diag(s1,s2,s3), s_i in {-1,+1}, so that
% diag(T*Q) is positive. The product T*M is preserved by setting Mpos = Q*M.
%
% Inputs:
%   T : 3x3 upper triangular
%   M : 3xN (optional, can be empty)
%
% Outputs:
%   Tpos : 3x3, T*Q with positive diagonal entries
%   Mpos : 3xN, Q*M
%   Q    : 3x3 diagonal reflection

d = diag(T);
s = sign(d);
s(s==0) = 1;
Q = diag(s);
Tpos = T * Q;

if nargin < 2 || isempty(M)
    Mpos = [];
else
    Mpos = Q * M;
end
end
