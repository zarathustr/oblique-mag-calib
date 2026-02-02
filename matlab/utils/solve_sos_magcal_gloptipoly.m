function out = solve_sos_magcal_gloptipoly(Yw, opts)
%SOLVE_SOS_MAGCAL_GLOPTIPOLY Globally optimize a short window via moment-SOS.
%
%   out = solve_sos_magcal_gloptipoly(Yw, opts)
%
% Solves (window length W = size(Yw,2)):
%   min_{T,h,{m_k}} sum_k ||Yw(:,k) - h - T m_k||_2^2
%   s.t. ||m_k||_2 = 1,  T upper triangular, diag(T) > 0
%
% This formulation matches the constrained ML problem on a window; the moment
% relaxation can be tight for small W and moderate noise, giving a global
% certificate.
%
% Requires GloptiPoly 3 + a supported SDP solver.
%
% Options:
%   opts.relax_order  (default [], let GloptiPoly choose)
%
% Outputs:
%   out.T, out.h, out.M
%   out.status, out.obj
%   out.moment_rank (if accessible)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'relax_order'), opts.relax_order = []; end

if exist('mpol', 'file') ~= 2
    error('GloptiPoly 3 not found on the MATLAB path.');
end

if exist('mset', 'file') == 2
    mset clear
end

W = size(Yw,2);

% Decision variables
mpol h 3 1
mpol a 3 1
mpol t12 t13 t23

eval(sprintf('mpol M 3 %d', W));

T = [a(1)^2, t12,    t13;
     0,      a(2)^2, t23;
     0,      0,      a(3)^2];

f = 0;
for k = 1:W
    rk = Yw(:,k) - h - T*M(:,k);
    f = f + rk.'*rk;
end

K = [];
for k = 1:W
    K = [K; M(:,k).'*M(:,k) == 1]; %#ok<AGROW>
end

if isempty(opts.relax_order)
    P = msdp(min(f), K);
else
    P = msdp(min(f), K, opts.relax_order);
end

[status, obj] = msol(P);

hsol = double(h);
Tsol = double(T);
Msol = double(M);

moment_rank = NaN;
try
    X = double(momentmatrix(P));
    s = svd(X);
    moment_rank = sum(s > 1e-6*s(1));
catch
    try
        X = double(mmat(P));
        s = svd(X);
        moment_rank = sum(s > 1e-6*s(1));
    catch
        X = [];
    end
end

out = struct();
out.T = Tsol;
out.h = hsol;
out.M = Msol;
out.status = status;
out.obj = obj;
out.moment_rank = moment_rank;
out.W = W;
end
