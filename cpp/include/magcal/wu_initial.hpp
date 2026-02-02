#pragma once

#include "magcal/types.hpp"

#include <Eigen/Eigenvalues>
#include <iostream>
#include <stdexcept>

namespace magcal {

struct WuInitResult {
    Mat3 A;   // ellipsoid matrix in centered form: (y-h)^T A (y-h) = 1
    Vec3 h;
    Mat3 R;   // upper triangular, R^T R = A
    Mat3 T;   // upper triangular, T = R^{-1}
};

inline WuInitResult wu_initial_estimate(const Mat3X& Y, bool verbose = true) {
    const int N = (int)Y.cols();
    if (N < 10) {
        throw std::runtime_error("Need at least 10 samples for ellipsoid fit");
    }

    Eigen::VectorXd x = Y.row(0).transpose();
    Eigen::VectorXd y = Y.row(1).transpose();
    Eigen::VectorXd z = Y.row(2).transpose();

    Eigen::MatrixXd D(N, 10);
    D.col(0) = x.array().square().matrix();
    D.col(1) = y.array().square().matrix();
    D.col(2) = z.array().square().matrix();
    D.col(3) = (2.0 * x.array() * y.array()).matrix();
    D.col(4) = (2.0 * x.array() * z.array()).matrix();
    D.col(5) = (2.0 * y.array() * z.array()).matrix();
    D.col(6) = x;
    D.col(7) = y;
    D.col(8) = z;
    D.col(9) = Eigen::VectorXd::Ones(N);

    Eigen::MatrixXd C = D.transpose() * D;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigen decomposition failed in ellipsoid fit");
    }

    // Smallest eigenvalue gives minimizer of ||D p|| with ||p||=1
    Eigen::VectorXd evals = es.eigenvalues();
    int idx = 0;
    evals.minCoeff(&idx);
    Eigen::VectorXd p = es.eigenvectors().col(idx);

    const double a11 = p(0);
    const double a22 = p(1);
    const double a33 = p(2);
    const double a12 = p(3);
    const double a13 = p(4);
    const double a23 = p(5);
    Vec3 b_e = p.segment<3>(6);
    const double c_e = p(9);

    Mat3 A_e;
    A_e << a11, a12, a13,
           a12, a22, a23,
           a13, a23, a33;
    A_e = 0.5 * (A_e + A_e.transpose());

    // Wu scaling: alpha = 4/(b_e^T A_e^{-1} b_e - 4 c_e)
    const double denom = (b_e.transpose() * A_e.inverse() * b_e)(0) - 4.0 * c_e;
    if (std::abs(denom) < 1e-16) {
        throw std::runtime_error("Degenerate ellipsoid fit: denom near zero");
    }
    const double alpha = 4.0 / denom;

    Mat3 A = alpha * A_e;
    Vec3 b = alpha * b_e;

    // Ensure SPD
    Eigen::SelfAdjointEigenSolver<Mat3> esA(A);
    if (esA.info() != Eigen::Success) {
        throw std::runtime_error("Eigen decomposition failed for A");
    }
    if (esA.eigenvalues().minCoeff() <= 0.0) {
        A = -A;
        b = -b;
    }

    Vec3 h = -0.5 * A.ldlt().solve(b);

    A = 0.5 * (A + A.transpose());
    Eigen::LLT<Mat3> llt(A);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("Cholesky factorization failed: A not SPD");
    }
    Mat3 R = llt.matrixU(); // upper triangular, R^T R = A

    WuInitResult out;
    out.A = A;
    out.h = h;
    out.R = R;
    out.T = R.inverse();

    if (verbose) {
        std::cout << "Wu init: N=" << N << "\n";
        std::cout << "h0=" << h.transpose() << "\n";
        std::cout << "T0=\n" << out.T << "\n";
    }
    return out;
}

} // namespace magcal
