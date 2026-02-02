#pragma once

#include "magcal/types.hpp"
#include "magcal/utils.hpp"
#include "magcal/wu_initial.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace magcal {

struct ManifoldOptions {
    int max_iters = 200;
    double grad_tol = 1e-7;
    double rel_cost_tol = 1e-12;
    double step0 = 1e-1;
    double backtrack_beta = 0.5;
    double armijo_c = 1e-4;
    int max_backtrack = 25;

    // Smoothness penalty on M (encourages continuous estimated directions)
    double lambda_smooth = 0.0;

    bool verbose = true;
};

struct ManifoldResult {
    Mat3 T;
    Vec3 h;
    Mat3X M;
    double cost = std::numeric_limits<double>::quiet_NaN();
    int iters = 0;
    bool converged = false;

    WuInitResult init;
};

struct State {
    Vec2 t2;   // unit
    Vec3 t3;   // unit
    Vec3 s;    // log-scales
    Vec3 h;
    Mat3X M;   // unit columns
};

inline Mat3 build_T(const State& x) {
    Vec3 u = x.s.array().exp().matrix();
    Mat3 T = Mat3::Zero();
    T(0,0) = u(0);
    T(0,1) = u(1) * x.t2(0);
    T(1,1) = u(1) * x.t2(1);
    T.col(2) = u(2) * x.t3;
    // Enforce strict upper triangular structure
    T(1,0) = 0.0; T(2,0) = 0.0; T(2,1) = 0.0;
    return T;
}

inline void T_to_params(const Mat3& T, Vec2& t2, Vec3& t3, Vec3& s) {
    const double u1 = T(0,0);
    const Eigen::Vector2d v2 = T.block<2,1>(0,1);
    double u2 = v2.norm();
    if (u2 < 1e-12) u2 = 1.0;
    t2 = v2 / u2;

    Eigen::Vector3d v3 = T.col(2);
    double u3 = v3.norm();
    if (u3 < 1e-12) u3 = 1.0;
    t3 = v3 / u3;

    s << std::log(std::max(1e-12, u1)), std::log(u2), std::log(u3);
}

inline double smoothness_cost(const Mat3X& M, double lambda) {
    if (lambda <= 0.0 || M.cols() < 2) return 0.0;
    double s = 0.0;
    for (int k = 0; k < M.cols()-1; ++k) {
        s += (M.col(k+1) - M.col(k)).squaredNorm();
    }
    return 0.5 * lambda * s;
}

inline void smoothness_grad(Mat3X& G, const Mat3X& M, double lambda) {
    if (lambda <= 0.0 || M.cols() < 2) return;
    const int N = M.cols();

    // k = 0
    G.col(0) += lambda * (M.col(0) - M.col(1));
    // k = N-1
    G.col(N-1) += lambda * (M.col(N-1) - M.col(N-2));
    // interior
    for (int k = 1; k < N-1; ++k) {
        G.col(k) += lambda * (2.0*M.col(k) - M.col(k-1) - M.col(k+1));
    }
}

inline double cost_function(const Mat3X& Y, const State& x, double lambda_smooth) {
    Mat3 T = build_T(x);
    const int N = Y.cols();

    Mat3X Res = Y - T * x.M - x.h * Eigen::RowVectorXd::Ones(N);
    double f = 0.5 * Res.squaredNorm();
    f += smoothness_cost(x.M, lambda_smooth);
    return f;
}

struct Grad {
    Vec2 gt2;
    Vec3 gt3;
    Vec3 gs;
    Vec3 gh;
    Mat3X gM;
};

inline Grad gradient_function(const Mat3X& Y, const State& x, double lambda_smooth) {
    const int N = Y.cols();
    Mat3 T = build_T(x);
    Mat3X Res = Y - T * x.M - x.h * Eigen::RowVectorXd::Ones(N);

    // Euclidean gradients wrt embedded variables
    Mat3 G_T = -Res * x.M.transpose();
    Vec3 G_h = -Res * Eigen::VectorXd::Ones(N);
    Mat3X G_M = -T.transpose() * Res;

    // Smoothness grad on M
    smoothness_grad(G_M, x.M, lambda_smooth);

    Vec3 u = x.s.array().exp().matrix();

    // Column 1: [u1;0;0]
    const double gs1 = u(0) * G_T(0,0);

    // Column 2: u2*[t2;0]
    const Eigen::Vector2d gcol2 = G_T.block<2,1>(0,1);
    const Vec2 gt2 = u(1) * gcol2;
    const double gs2 = u(1) * (x.t2.dot(gcol2));

    // Column 3: u3*t3
    const Vec3 gt3 = u(2) * G_T.col(2);
    const double gs3 = u(2) * (x.t3.dot(G_T.col(2)));

    Grad g;
    g.gt2 = gt2;
    g.gt3 = gt3;
    g.gs  << gs1, gs2, gs3;
    g.gh  = G_h;
    g.gM  = std::move(G_M);
    return g;
}

inline Vec2 proj_sphere(const Vec2& t, const Vec2& g) {
    return g - (t.dot(g)) * t;
}
inline Vec3 proj_sphere(const Vec3& t, const Vec3& g) {
    return g - (t.dot(g)) * t;
}

inline void proj_oblique_inplace(const Mat3X& M, Mat3X& G) {
    // project each column to tangent space of sphere
    for (int k = 0; k < M.cols(); ++k) {
        const double dot = M.col(k).dot(G.col(k));
        G.col(k) -= dot * M.col(k);
    }
}

inline double grad_norm(const Vec2& gt2, const Vec3& gt3, const Vec3& gs, const Vec3& gh, const Mat3X& gM) {
    double n2 = gt2.squaredNorm() + gt3.squaredNorm() + gs.squaredNorm() + gh.squaredNorm();
    n2 += gM.squaredNorm();
    return std::sqrt(n2);
}

inline State retract(const State& x, const Vec2& dt2, const Vec3& dt3, const Vec3& ds, const Vec3& dh, const Mat3X& dM, double step) {
    State y;
    y.s = x.s + step * ds;
    y.h = x.h + step * dh;

    y.t2 = x.t2 + step * dt2;
    if (y.t2.norm() < 1e-12) y.t2 = Vec2(1.0, 0.0);
    y.t2.normalize();

    y.t3 = x.t3 + step * dt3;
    if (y.t3.norm() < 1e-12) y.t3 = Vec3(1.0, 0.0, 0.0);
    y.t3.normalize();

    y.M = x.M + step * dM;
    y.M = normalize_columns(y.M);
    return y;
}

inline ManifoldResult solve_product_manifold(const Mat3X& Y, const ManifoldOptions& opt = {}) {
    const int N = Y.cols();

    // 1) Wu initialization
    WuInitResult init = wu_initial_estimate(Y, opt.verbose);

    // 2) Initial state
    State x;
    T_to_params(init.T, x.t2, x.t3, x.s);
    x.h = init.h;

    // initial M from affine inverse
    Mat3X M0 = init.T.inverse() * (Y - x.h * Eigen::RowVectorXd::Ones(N));
    x.M = normalize_columns(M0);

    // Ensure unit norms
    if (x.t2.norm() < 1e-12) x.t2 = Vec2(1.0, 0.0);
    x.t2.normalize();
    if (x.t3.norm() < 1e-12) x.t3 = Vec3(1.0, 0.0, 0.0);
    x.t3.normalize();

    double f = cost_function(Y, x, opt.lambda_smooth);

    if (opt.verbose) {
        std::cout << "Initial cost: " << f << "\n";
    }

    double step0 = opt.step0;

    for (int iter = 0; iter < opt.max_iters; ++iter) {
        Grad g = gradient_function(Y, x, opt.lambda_smooth);

        // Project to tangent spaces
        Vec2 gt2_tan = proj_sphere(x.t2, g.gt2);
        Vec3 gt3_tan = proj_sphere(x.t3, g.gt3);
        Mat3X gM_tan = g.gM;
        proj_oblique_inplace(x.M, gM_tan);

        const double gn = grad_norm(gt2_tan, gt3_tan, g.gs, g.gh, gM_tan);
        if (opt.verbose && (iter % 10 == 0 || iter == opt.max_iters-1)) {
            std::cout << "iter " << iter << ": cost=" << f << " gradnorm=" << gn << "\n";
        }

        if (gn < opt.grad_tol) {
            ManifoldResult res;
            res.T = build_T(x);
            res.h = x.h;
            res.M = x.M;
            res.cost = f;
            res.iters = iter;
            res.converged = true;
            res.init = init;
            return res;
        }

        // Descent direction is negative projected gradient
        Vec2 dt2 = -gt2_tan;
        Vec3 dt3 = -gt3_tan;
        Vec3 ds  = -g.gs;
        Vec3 dh  = -g.gh;
        Mat3X dM = -gM_tan;

        // Backtracking Armijo line search
        double step = step0;
        const double gdot = gn*gn; // approx ||grad||^2 in Riemannian metric
        double f_new = std::numeric_limits<double>::infinity();
        State x_new;

        for (int bt = 0; bt < opt.max_backtrack; ++bt) {
            x_new = retract(x, dt2, dt3, ds, dh, dM, step);
            f_new = cost_function(Y, x_new, opt.lambda_smooth);
            if (f_new <= f - opt.armijo_c * step * gdot) {
                break;
            }
            step *= opt.backtrack_beta;
        }

        const double rel_drop = (f - f_new) / std::max(1.0, std::abs(f));

        x = std::move(x_new);
        f = f_new;

        // Simple step size adaptation
        if (rel_drop > 1e-3) {
            step0 = std::min(1.0, step * 1.2);
        } else {
            step0 = std::max(1e-6, step * 0.8);
        }

        if (rel_drop < opt.rel_cost_tol) {
            ManifoldResult res;
            res.T = build_T(x);
            res.h = x.h;
            res.M = x.M;
            res.cost = f;
            res.iters = iter + 1;
            res.converged = true;
            res.init = init;
            return res;
        }
    }

    ManifoldResult res;
    res.T = build_T(x);
    res.h = x.h;
    res.M = x.M;
    res.cost = f;
    res.iters = opt.max_iters;
    res.converged = false;
    res.init = init;
    return res;
}

} // namespace magcal
