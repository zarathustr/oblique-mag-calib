#pragma once

#include "magcal/types.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace magcal {

inline Mat3 ellipsoid_A_from_T(const Mat3& T) {
    Mat3 Ti = T.inverse();
    return Ti.transpose() * Ti; // SPD
}

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Eigen::Vector3i> faces; // 0-indexed triangles
};

inline Mesh make_ellipsoid_mesh(const Mat3& A, const Vec3& h, int nu = 60, int nv = 30) {
    Eigen::SelfAdjointEigenSolver<Mat3> es(A);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigen decomposition failed in make_ellipsoid_mesh");
    }
    const Vec3 evals = es.eigenvalues();
    const Mat3 Q = es.eigenvectors();
    if ((evals.array() <= 0.0).any()) {
        throw std::runtime_error("Ellipsoid matrix A must be SPD (positive eigenvalues)");
    }

    Vec3 axes;
    axes(0) = 1.0 / std::sqrt(evals(0));
    axes(1) = 1.0 / std::sqrt(evals(1));
    axes(2) = 1.0 / std::sqrt(evals(2));

    const int Nu = std::max(3, nu);
    const int Nv = std::max(2, nv);

    Mesh mesh;
    mesh.vertices.reserve((Nu+1) * (Nv+1));

    const double pi = 3.141592653589793238462643383279502884;

    for (int iv = 0; iv <= Nv; ++iv) {
        const double v = pi * (double)iv / (double)Nv; // 0..pi
        const double sv = std::sin(v);
        const double cv = std::cos(v);
        for (int iu = 0; iu <= Nu; ++iu) {
            const double u = 2.0 * pi * (double)iu / (double)Nu; // 0..2pi
            const double cu = std::cos(u);
            const double su = std::sin(u);

            Vec3 p0;
            p0(0) = axes(0) * cu * sv;
            p0(1) = axes(1) * su * sv;
            p0(2) = axes(2) * cv;

            Vec3 p = Q * p0 + h;
            mesh.vertices.push_back(p);
        }
    }

    // faces: connect grid cells
    auto idx = [Nu](int iu, int iv) {
        return iv * (Nu+1) + iu;
    };

    for (int iv = 0; iv < Nv; ++iv) {
        for (int iu = 0; iu < Nu; ++iu) {
            int i00 = idx(iu, iv);
            int i10 = idx(iu+1, iv);
            int i01 = idx(iu, iv+1);
            int i11 = idx(iu+1, iv+1);

            // two triangles
            mesh.faces.push_back(Eigen::Vector3i(i00, i10, i11));
            mesh.faces.push_back(Eigen::Vector3i(i00, i11, i01));
        }
    }

    return mesh;
}

inline void write_obj(const Mesh& mesh, const std::string& path) {
    std::ofstream fout(path);
    if (!fout.good()) {
        throw std::runtime_error("Failed to open for writing: " + path);
    }

    for (const auto& v : mesh.vertices) {
        fout << "v " << v(0) << " " << v(1) << " " << v(2) << "\n";
    }
    for (const auto& f : mesh.faces) {
        // OBJ is 1-indexed
        fout << "f " << (f(0)+1) << " " << (f(1)+1) << " " << (f(2)+1) << "\n";
    }
}

inline void write_obj_ellipsoid(const Mat3& A, const Vec3& h, const std::string& path, int nu = 60, int nv = 30) {
    Mesh mesh = make_ellipsoid_mesh(A, h, nu, nv);
    write_obj(mesh, path);
}

} // namespace magcal
