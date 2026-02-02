#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

namespace magcal {

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

inline bool contains_case_insensitive(const std::string& haystack, const std::string& needle) {
    return to_lower(haystack).find(to_lower(needle)) != std::string::npos;
}

inline std::vector<std::string> split_csv(const std::string& line, char sep = ',') {
    std::vector<std::string> out;
    std::string token;
    std::stringstream ss(line);
    while (std::getline(ss, token, sep)) {
        out.push_back(token);
    }
    // If line ends with a separator, getline will miss the last empty token.
    if (!line.empty() && line.back() == sep) {
        out.emplace_back("");
    }
    return out;
}

inline void normalize_columns_inplace(Eigen::Ref<Eigen::MatrixXd> X, double eps = 1e-12) {
    for (int k = 0; k < X.cols(); ++k) {
        double n = X.col(k).norm();
        if (n < eps) n = 1.0;
        X.col(k) /= n;
    }
}

inline Eigen::Matrix<double,3,Eigen::Dynamic> normalize_columns(const Eigen::Matrix<double,3,Eigen::Dynamic>& X, double eps = 1e-12) {
    Eigen::Matrix<double,3,Eigen::Dynamic> Y = X;
    for (int k = 0; k < Y.cols(); ++k) {
        double n = Y.col(k).norm();
        if (n < eps) n = 1.0;
        Y.col(k) /= n;
    }
    return Y;
}

} // namespace magcal
