#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace magcal {

using Mat3X = Eigen::Matrix<double, 3, Eigen::Dynamic>;
using Vec3  = Eigen::Matrix<double, 3, 1>;
using Vec2  = Eigen::Matrix<double, 2, 1>;
using Mat3  = Eigen::Matrix<double, 3, 3>;

struct MagData {
    Mat3X Y;                 // 3xN magnetometer measurements
    std::vector<double> t;   // optional timestamps (seconds), can be empty
    std::vector<std::string> source_files;
};

} // namespace magcal
