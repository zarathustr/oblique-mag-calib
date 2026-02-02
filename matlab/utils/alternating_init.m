function [T,h,M,info] = alternating_init(Y, maxiter)
%ALTERNATING_INIT Alternating minimization using closed-form block updates.
% Updates:
% 1) M_k = normalize(T^{-1}(y_k - h))
% 2) (T,h) from linear least squares given M

if nargin < 2, maxiter = 50; end
N = size(Y,2);

h = mean(Y,2);
T = eye(3);

prev = inf;
for it = 1:maxiter
    A = T \ (Y - h*ones(1,N));
    M = normalize_columns(A);
    [T,h] = solve_ls_T_h_given_M(Y,M);

    E = Y - h*ones(1,N) - T*M;
    cost = 0.5*sum(E(:).^2)/N;
    if abs(prev - cost) < 1e-12*max(1,prev), break; end
    prev = cost;
end

info.iters = it;
info.cost = cost;
end
