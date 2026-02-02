function [T,h,M,info] = alternating_refine_from_seed(Y, T, h, maxiter)
%ALTERNATING_REFINE_FROM_SEED Alternating minimization from a provided seed.
% Updates:
%   1) M_k = normalize(T^{-1}(y_k - h))
%   2) (T,h) from linear least squares given M
%
% Inputs:
%   Y 3xN samples
%   T 3x3 (should be upper triangular)
%   h 3x1
%   maxiter (default 50)

if nargin < 4 || isempty(maxiter)
    maxiter = 50;
end
N = size(Y,2);
oneN = ones(1,N);

prev = inf;
for it = 1:maxiter
    A = T \ (Y - h*oneN);
    M = normalize_columns(A);
    [T, h] = solve_ls_T_h_given_M(Y, M);

    E = Y - h*oneN - T*M;
    cost = 0.5*sum(E(:).^2)/N;
    if abs(prev - cost) < 1e-12*max(1,prev)
        break;
    end
    prev = cost;
end

info.iters = it;
info.cost = cost;
end
