#pragma once

#include "magcal/types.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace magcal {

// Loader for UST Nautilus IMU logs stored as plain text.
//
// Expected format: each non-empty line contains 10 numeric fields
// separated by whitespace (tabs or spaces):
//
//   seq  gx gy gz  ax ay az  mx my mz
//
// The magnetometer measurement is taken as the last three fields.
// The sequence number can optionally be stored in MagData::t.

struct Txt10LoadOptions {
    int downsample = 1;          // keep every downsample-th row
    bool verbose = true;
    bool store_seq_in_t = true;  // store the seq index as t[k]
};

inline MagData load_nautilus_txt10(const std::string& filepath, const Txt10LoadOptions& opt = {}) {
    std::ifstream fin(filepath);
    if (!fin.good()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }

    std::string line;
    std::vector<Vec3> samples;
    std::vector<double> t;
    int row_idx = 0;

    while (std::getline(fin, line)) {
        // Trim quickly
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        long long seq = 0;
        double gx = 0, gy = 0, gz = 0;
        double ax = 0, ay = 0, az = 0;
        double mx = 0, my = 0, mz = 0;

        if (!(ss >> seq >> gx >> gy >> gz >> ax >> ay >> az >> mx >> my >> mz)) {
            // Skip malformed lines
            continue;
        }

        if (opt.downsample <= 1 || (row_idx % opt.downsample) == 0) {
            Vec3 v;
            v << mx, my, mz;
            samples.push_back(v);
            if (opt.store_seq_in_t) {
                t.push_back(static_cast<double>(seq));
            }
        }
        row_idx++;
    }

    if (samples.empty()) {
        throw std::runtime_error("No magnetometer samples loaded from: " + filepath);
    }

    Mat3X Y(3, static_cast<int>(samples.size()));
    for (int k = 0; k < static_cast<int>(samples.size()); ++k) {
        Y.col(k) = samples[static_cast<size_t>(k)];
    }

    MagData out;
    out.Y = std::move(Y);
    if (opt.store_seq_in_t) {
        out.t = std::move(t);
    }
    out.source_files = {filepath};

    if (opt.verbose) {
        std::cout << "Loaded " << out.Y.cols() << " mag samples from " << filepath
                  << " (txt10, downsample=" << opt.downsample << ")\n";
    }
    return out;
}

// Convenience: load a directory containing a mixture of .csv and .txt logs.
// - .csv files are parsed with load_magnetometer_csv (Microstrain-style)
// - .txt files are parsed with load_nautilus_txt10
template <typename CsvLoaderFn>
inline MagData load_magnetometer_directory_mixed(
    const std::string& dirpath,
    CsvLoaderFn csv_loader,
    const Txt10LoadOptions& txt_opt = {})
{
    namespace fs = std::filesystem;
    if (!fs::exists(dirpath)) {
        throw std::runtime_error("Directory does not exist: " + dirpath);
    }

    std::vector<fs::path> files;
    for (auto& p : fs::directory_iterator(dirpath)) {
        if (!p.is_regular_file()) continue;
        auto ext = p.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (ext == ".csv" || ext == ".txt") {
            files.push_back(p.path());
        }
    }

    if (files.empty()) {
        throw std::runtime_error("No .csv or .txt files found in: " + dirpath);
    }

    std::sort(files.begin(), files.end());

    std::vector<Vec3> all;
    std::vector<double> t_all;
    std::vector<std::string> src;

    for (const auto& p : files) {
        auto ext = p.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (ext == ".csv") {
            MagData d = csv_loader(p.string());
            src.push_back(p.string());
            all.reserve(all.size() + static_cast<size_t>(d.Y.cols()));
            for (int k = 0; k < d.Y.cols(); ++k) {
                all.push_back(d.Y.col(k));
            }
            // ignore timestamps for CSV, as they are not standardized
        } else if (ext == ".txt") {
            MagData d = load_nautilus_txt10(p.string(), txt_opt);
            src.push_back(p.string());
            all.reserve(all.size() + static_cast<size_t>(d.Y.cols()));
            for (int k = 0; k < d.Y.cols(); ++k) {
                all.push_back(d.Y.col(k));
            }
            if (!d.t.empty()) {
                t_all.reserve(t_all.size() + d.t.size());
                t_all.insert(t_all.end(), d.t.begin(), d.t.end());
            }
        }
    }

    Mat3X Y(3, static_cast<int>(all.size()));
    for (int k = 0; k < static_cast<int>(all.size()); ++k) {
        Y.col(k) = all[static_cast<size_t>(k)];
    }

    MagData out;
    out.Y = std::move(Y);
    out.source_files = std::move(src);
    if (!t_all.empty()) {
        out.t = std::move(t_all);
    }
    return out;
}

} // namespace magcal
