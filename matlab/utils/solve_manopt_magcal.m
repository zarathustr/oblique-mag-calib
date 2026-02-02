function out = solve_manopt_magcal(Y, opts)
%SOLVE_MANOPT_MAGCAL Calibrate magnetometer using Manopt on a product manifold.
%
%   out = solve_manopt_magcal(Y, opts)
%
% Measurement model:
%   Y(:,k) = h + T * M(:,k) + noise
% with ||M(:,k)||_2 = 1 (columns on oblique manifold) and T upper triangular.
%
% This uses the same product-manifold parameterization as demo_manopt_product_manifold:
%   t11           = exp(sigma1)
%   [t12;t22]     = exp(sigma2) * u2, u2 in S^1
%   [t13;t23;t33] = exp(sigma3) * u3, u3 in S^2
%
% Optional smoothness penalty encourages continuous motion:
%   (lambda/2) * sum_k ||M(:,k+1) - M(:,k)||^2
%
% Inputs:
%   Y    3xN raw magnetometer samples
%   opts struct with optional fields:
%        lambda          smoothness weight (default 1e-2)
%        maxiter         trust-regions iterations (default 200)
%        tolgradnorm     TR stopping tolerance (default 1e-8)
%        init_maxiter    alternating initialization iterations (default 50)
%        verbosity       Manopt verbosity (default 2)
%
% Output struct:
%   out.T, out.h, out.M     estimated parameters
%   out.T0, out.h0, out.M0  initialization
%   out.info0               init diagnostics
%   out.info                Manopt info array
%   out.cost_final          final cost
%   out.options             Manopt options used

if nargin < 2, opts = struct(); end
if ~isfield(opts,'lambda'), opts.lambda = 1e-2; end
if ~isfield(opts,'maxiter'), opts.maxiter = 200; end
if ~isfield(opts,'tolgradnorm'), opts.tolgradnorm = 1e-8; end
if ~isfield(opts,'init_maxiter'), opts.init_maxiter = 50; end
if ~isfield(opts,'verbosity'), opts.verbosity = 2; end
if ~isfield(opts,'init_method'), opts.init_method = 'wu_alt'; end

N = size(Y,2);
oneN = ones(1,N);

% Initialization
switch lower(opts.init_method)
    case {'alternating','alt'}
        [T0, h0, M0, info0] = alternating_init(Y, opts.init_maxiter);
    case {'wu','wu_alt','wu+alt'}
        [Tinit, hinit] = wu2015_initial_T_h(Y);
        [T0, h0, M0, info0] = alternating_refine_from_seed(Y, Tinit, hinit, opts.init_maxiter);
    otherwise
        error('Unknown opts.init_method: %s', opts.init_method);
end
[sigma0, u20, u30] = triU_to_oblique_params(T0);

% Product manifold
manifold = productmanifold(struct( ...
    'sigma', euclideanfactory(3,1), ...
    'u2', spherefactory(2), ...
    'u3', spherefactory(3), ...
    'h', euclideanfactory(3,1), ...
    'M', obliquefactory(3,N)));

problem.M = manifold;

lambda = opts.lambda;
problem.cost = @cost;
problem.egrad = @egrad;

x0 = struct('sigma', sigma0, 'u2', u20, 'u3', u30, 'h', h0, 'M', M0);

options = struct();
options.verbosity = opts.verbosity;
options.maxiter = opts.maxiter;
options.tolgradnorm = opts.tolgradnorm;

[xopt, cost_final, info, options] = trustregions(problem, x0, options);

T = build_T_from_vars(xopt.sigma, xopt.u2, xopt.u3);
[T, M, Q] = enforce_positive_diag_gauge(T, xopt.M);
h = xopt.h;

out = struct();
out.T = T;
out.h = h;
out.M = M;
out.T0 = T0;
out.h0 = h0;
out.M0 = M0;
out.info0 = info0;
out.info = info;
out.cost_final = cost_final;
out.options = options;
out.lambda = lambda;

    function f = cost(x)
        Tx = build_T_from_vars(x.sigma, x.u2, x.u3);
        E = Y - x.h*oneN - Tx*x.M;
        f = 0.5*sum(E(:).^2);
        if lambda > 0
            dM = diff(x.M,1,2);
            f = f + 0.5*lambda*sum(dM(:).^2);
        end
    end

    function g = egrad(x)
        Tx = build_T_from_vars(x.sigma, x.u2, x.u3);
        E = Y - x.h*oneN - Tx*x.M;

        gT = -E*x.M.';    % d/dT
        gh = -E*oneN.';   % d/dh
        gM = -(Tx.')*E;   % d/dM

        if lambda > 0
            gM(:,1) = gM(:,1) + lambda*(x.M(:,1) - x.M(:,2));
            gM(:,N) = gM(:,N) + lambda*(x.M(:,N) - x.M(:,N-1));
            if N > 2
                gM(:,2:N-1) = gM(:,2:N-1) + lambda*(2*x.M(:,2:N-1) - x.M(:,1:N-2) - x.M(:,3:N));
            end
        end

        % Chain rule from gT to (sigma,u2,u3)
        t11 = exp(x.sigma(1));
        gs1 = gT(1,1) * t11;

        u2 = x.u2(:);
        a2 = exp(x.sigma(2))*u2;
        gcol2 = [gT(1,2); gT(2,2)];
        gs2 = dot(gcol2, a2);
        gu2 = exp(x.sigma(2)) * gcol2;

        u3 = x.u3(:);
        a3 = exp(x.sigma(3))*u3;
        gcol3 = [gT(1,3); gT(2,3); gT(3,3)];
        gs3 = dot(gcol3, a3);
        gu3 = exp(x.sigma(3)) * gcol3;

        g = struct();
        g.sigma = [gs1; gs2; gs3];
        g.u2 = gu2;
        g.u3 = gu3;
        g.h = gh;
        g.M = gM;
    end
end
