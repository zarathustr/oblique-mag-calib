%% run_nautilus_experiment.m
% Real-data experiment using UST Nautilus TXT logs.
%
% The TXT format is:
%   seq gx gy gz ax ay az mx my mz
%
% This script:
%   1) loads and concatenates magnetometer samples from ../data/nautilus/*.txt
%   2) runs the Manopt product-manifold solver (oblique product manifold)
%   3) runs a windowed global moment-SOS solver (GloptiPoly) if available
%   4) plots fitted ellipsoids in raw measurement space using ellipsoid()
%      and rotating the mesh to match the fitted quadratic form

clear; close all; clc;

thisdir = fileparts(mfilename('fullpath'));
root = fullfile(thisdir, '..');
addpath(genpath(root));

% -------------------------
% 1) Load data
% -------------------------
dataFolder = fullfile(root, 'data', 'nautilus');
data = load_nautilus_dataset(dataFolder);
Yfull = data.Y;
fprintf('Loaded %d samples from %d files\n', data.N, numel(data.files));

% Downsample for speed (keeps time-series continuity)
ds = 2; % set ds=1 for full resolution
Y = Yfull(:, 1:ds:end);
N = size(Y,2);
oneN = ones(1,N);
fprintf('Downsample factor = %d, working N = %d\n', ds, N);

% -------------------------
% 2) Initialization ellipsoid (Wu 2015)
% -------------------------
[Twu, hwu] = wu2015_initial_T_h(Y);
Awu = A_from_T(Twu);

% -------------------------
% 3) Manopt product-manifold solver
% -------------------------
doManopt = exist('trustregions','file') == 2;
if ~doManopt
    error('Manopt trustregions() not found. Add Manopt to your MATLAB path.');
end

optsM = struct();
optsM.lambda = 1e-2;
optsM.maxiter = 200;
optsM.tolgradnorm = 1e-8;
optsM.init_maxiter = 50;
optsM.init_method = 'wu_alt';
optsM.verbosity = 2;

outM = solve_manopt_magcal(Y, optsM);
Aman = A_from_T(outM.T);
hman = outM.h;

% -------------------------
% 4) Global SOS solver on short windows
% -------------------------
doSOS = exist('mpol','file') == 2;
Tsos = []; hsos = []; Asos = [];
SOS_details = struct();

if doSOS
    W = 12;           % window length (increase cautiously)
    numWins = 10;     % number of windows to solve
    starts = round(linspace(1, N-W+1, numWins));

    A_list = zeros(3,3,numWins);
    h_list = zeros(3,numWins);
    status_list = zeros(1,numWins);
    rank_list = nan(1,numWins);
    obj_list = nan(1,numWins);

    fprintf('\nRunning SOS on %d windows, W=%d...\n', numWins, W);
    for i = 1:numWins
        idx = starts(i);
        Yw = Y(:, idx:idx+W-1);
        try
            outS = solve_sos_magcal_gloptipoly(Yw);
            status_list(i) = outS.status;
            rank_list(i) = outS.moment_rank;
            obj_list(i) = outS.obj;
            A_list(:,:,i) = A_from_T(outS.T);
            h_list(:,i) = outS.h;
            fprintf('  window %2d/%2d start=%5d: status=%d  obj=%.3e  rank=%g\n', ...
                i, numWins, idx, outS.status, outS.obj, outS.moment_rank);
        catch ME
            warning('SOS window %d failed: %s', i, ME.message);
            status_list(i) = -999;
        end
    end

    good = (status_list == 1) | (status_list == 0);
    if any(good)
        Asos = mean(A_list(:,:,good), 3);
        hsos = mean(h_list(:,good), 2);
        % Reconstruct an upper-triangular T from A = T^{-T}T^{-1}.
        U = chol(Asos, 'upper');
        Tsos = inv(U);
    else
        warning('No successful SOS windows. Skipping SOS ellipsoid.');
    end

    SOS_details.starts = starts;
    SOS_details.status = status_list;
    SOS_details.rank = rank_list;
    SOS_details.obj = obj_list;
else
    fprintf('\nGloptiPoly not found on path; skipping SOS solver.\n');
end

% -------------------------
% 5) Plot raw cloud + fitted ellipsoids
% -------------------------
figure('Name','Nautilus raw data with fitted ellipsoids');
plot3(Y(1, :), Y(2, :), Y(3, :), '.', 'MarkerSize', 16); hold on
axis equal
grid minor
xlabel('$m_x$', 'FontSize', 18, 'Interpreter', 'LaTeX'); 
ylabel('$m_y$', 'FontSize', 18, 'Interpreter', 'LaTeX'); 
zlabel('$m_z$', 'FontSize', 18, 'Interpreter', 'LaTeX');
title('Raw Data and Fitted Ellipsoids', 'FontSize', 18, 'Interpreter', 'LaTeX');

% Wu init ellipsoid (uses ellipsoid() + rotation)
[Xw, Yw_, Zw] = ellipsoid_mesh_from_Ah(Awu, hwu, 50);
sw = surf(Xw, Yw_, Zw);
set(sw, 'FaceAlpha', 1.0, 'EdgeAlpha', 1);

% Manopt ellipsoid
[Xm, Ym, Zm] = ellipsoid_mesh_from_Ah(Aman, hman, 60);
sm = surf(Xm, Ym, Zm);
set(sm, 'FaceAlpha', 0.03, 'EdgeAlpha', 1, 'EdgeColor', "#A2142F", 'FaceColor', "#A2142F");

lg = {'Raw Samples', 'Wu init ellipsoid','Manopt ellipsoid'};

% SOS ellipsoid (if available)
if ~isempty(Asos)
    [Xs, Ys, Zs] = ellipsoid_mesh_from_Ah(Asos, hsos, 60);
    ss = surf(Xs, Ys, Zs);
    set(ss, 'FaceAlpha', 0.10, 'EdgeAlpha', 0.10);
    lg{end+1} = 'SOS ellipsoid (avg windows)';
end
legend(lg, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'LaTeX');

% -------------------------
% 6) Plot calibrated points (should be close to a sphere)
% -------------------------
Mcal = outM.T \ (Y - outM.h*oneN);
Mnorm = vecnorm(Mcal, 2, 1);

figure('Name', 'Nautilus calibrated data');
subplot(1, 2, 1);
plot(Mnorm); 
grid minor
xlabel('Sample Index', 'FontSize', 18, 'Interpreter', 'LaTeX'); 
ylabel('$||\mathbf{m}_k||$', 'FontSize', 18, 'Interpreter', 'LaTeX');
ylim([0.9, 1.1])
title('\texttt{Manopt} Calibrated Norms', 'FontSize', 18, 'Interpreter', 'LaTeX');

subplot(1, 2, 2);
plot3(Mcal(1, :), Mcal(2, :), Mcal(3, :), '.', ...
    'MarkerSize', 16, 'Color', [0.7, 0.2, 0.2]); hold on
axis equal
grid minor
xlabel('$m_x$', 'FontSize', 18, 'Interpreter', 'LaTeX'); 
ylabel('$m_y$', 'FontSize', 18, 'Interpreter', 'LaTeX'); 
zlabel('$m_z$', 'FontSize', 18, 'Interpreter', 'LaTeX');
title('Calibrated Magnetometer Point Cloud', 'FontSize', 18, 'Interpreter', 'LaTeX');

% Draw unit sphere for reference
[XS, YS, ZS] = sphere(30);
surf(XS, YS, ZS, 'FaceAlpha', 0.6, 'EdgeAlpha', 1);

% -------------------------
% 7) Simple quantitative summaries
% -------------------------
q_wu = sum((Y - hwu*oneN) .* (Awu*(Y - hwu*oneN)), 1);
q_man = sum((Y - hman*oneN) .* (Aman*(Y - hman*oneN)), 1);

fprintf('\nQuadratic-form statistics (target 1):\n');
fprintf('  Wu init:   mean=%.4f  std=%.4f\n', mean(q_wu), std(q_wu));
fprintf('  Manopt:    mean=%.4f  std=%.4f\n', mean(q_man), std(q_man));
if ~isempty(Asos)
    q_sos = sum((Y - hsos*oneN) .* (Asos*(Y - hsos*oneN)), 1);
    fprintf('  SOS (avg): mean=%.4f  std=%.4f\n', mean(q_sos), std(q_sos));
end

fprintf('\nDone.\n');
