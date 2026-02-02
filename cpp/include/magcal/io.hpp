#pragma once

#include "magcal/types.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace magcal {

inline void write_calibrated_csv(const Mat3X& Mcal, const std::string& path) {
    std::ofstream f(path);
    if (!f.good()) {
        throw std::runtime_error("Failed to open for writing: " + path);
    }
    f << "mx,my,mz\n";
    f.setf(std::ios::fixed);
    f << std::setprecision(10);
    for (int k = 0; k < Mcal.cols(); ++k) {
        f << Mcal(0,k) << "," << Mcal(1,k) << "," << Mcal(2,k) << "\n";
    }
}

template <typename Derived>
inline void write_vector_txt(const Eigen::MatrixBase<Derived>& v, const std::string& path) {
    std::ofstream f(path);
    if (!f.good()) {
        throw std::runtime_error("Failed to open for writing: " + path);
    }
    f.setf(std::ios::fixed);
    f << std::setprecision(12);
    for (int i = 0; i < v.size(); ++i) {
        f << v.derived()(i) << "\n";
    }
}

template <typename Derived>
inline void write_matrix_txt(const Eigen::MatrixBase<Derived>& M, const std::string& path) {
    std::ofstream f(path);
    if (!f.good()) {
        throw std::runtime_error("Failed to open for writing: " + path);
    }
    f.setf(std::ios::fixed);
    f << std::setprecision(12);
    for (int r = 0; r < M.rows(); ++r) {
        for (int c = 0; c < M.cols(); ++c) {
            f << M.derived()(r,c);
            if (c + 1 < M.cols()) f << " ";
        }
        f << "\n";
    }
}

} // namespace magcal
