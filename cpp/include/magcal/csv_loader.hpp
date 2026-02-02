#pragma once

#include "magcal/types.hpp"
#include "magcal/utils.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace magcal {

struct CsvLoadOptions {
    int downsample = 1;                 // keep every downsample-th row
    std::string x_name = "x mag";      // case-insensitive substring match
    std::string y_name = "y mag";
    std::string z_name = "z mag";
    bool verbose = true;
};

inline MagData load_magnetometer_csv(const std::string& filepath, const CsvLoadOptions& opt = {}) {
    std::ifstream fin(filepath);
    if (!fin.good()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }

    std::string line;
    bool in_data = false;
    bool header_read = false;
    int ix = -1, iy = -1, iz = -1;

    std::vector<Vec3> samples;
    int row_idx = 0;

    while (std::getline(fin, line)) {
        if (!in_data) {
            if (contains_case_insensitive(line, "DATA_START")) {
                in_data = true;
                continue;
            }
            continue;
        }

        if (in_data && !header_read) {
            // This line is the CSV header
            auto hdr = split_csv(line, ',');
            for (int i = 0; i < (int)hdr.size(); ++i) {
                const std::string h = to_lower(hdr[i]);
                if (ix < 0 && h.find(to_lower(opt.x_name)) != std::string::npos) ix = i;
                if (iy < 0 && h.find(to_lower(opt.y_name)) != std::string::npos) iy = i;
                if (iz < 0 && h.find(to_lower(opt.z_name)) != std::string::npos) iz = i;
            }
            if (ix < 0 || iy < 0 || iz < 0) {
                std::stringstream ss;
                ss << "Could not locate magnetometer columns in header for file: " << filepath
                   << "\nHeader was: " << line;
                throw std::runtime_error(ss.str());
            }
            header_read = true;
            continue;
        }

        if (line.empty()) continue;
        auto tok = split_csv(line, ',');
        if (ix >= (int)tok.size() || iy >= (int)tok.size() || iz >= (int)tok.size()) continue;

        if (opt.downsample <= 1 || (row_idx % opt.downsample) == 0) {
            try {
                // Some files have leading empty GPS fields; stod handles leading spaces but not empty strings.
                if (tok[ix].empty() || tok[iy].empty() || tok[iz].empty()) {
                    // skip rows with missing mag
                } else {
                    Vec3 v;
                    v(0) = std::stod(tok[ix]);
                    v(1) = std::stod(tok[iy]);
                    v(2) = std::stod(tok[iz]);
                    samples.push_back(v);
                }
            } catch (...) {
                // Skip malformed line
            }
        }
        row_idx++;
    }

    if (samples.empty()) {
        throw std::runtime_error("No magnetometer samples loaded from: " + filepath);
    }

    Mat3X Y(3, (int)samples.size());
    for (int k = 0; k < (int)samples.size(); ++k) {
        Y.col(k) = samples[k];
    }

    MagData out;
    out.Y = std::move(Y);
    out.source_files = {filepath};
    if (opt.verbose) {
        std::cout << "Loaded " << out.Y.cols() << " mag samples from " << filepath
                  << " (downsample=" << opt.downsample << ")\n";
    }
    return out;
}

inline MagData load_magnetometer_directory(const std::string& dirpath, const CsvLoadOptions& opt = {}) {
    namespace fs = std::filesystem;
    if (!fs::exists(dirpath)) {
        throw std::runtime_error("Directory does not exist: " + dirpath);
    }

    std::vector<fs::path> csvs;
    for (auto& p : fs::directory_iterator(dirpath)) {
        if (!p.is_regular_file()) continue;
        const auto ext = to_lower(p.path().extension().string());
        if (ext == ".csv") {
            csvs.push_back(p.path());
        }
    }

    if (csvs.empty()) {
        throw std::runtime_error("No .csv files found in: " + dirpath);
    }

    std::sort(csvs.begin(), csvs.end());

    std::vector<Vec3> all;
    std::vector<std::string> files;

    for (const auto& p : csvs) {
        auto d = load_magnetometer_csv(p.string(), opt);
        files.push_back(p.string());
        all.reserve(all.size() + d.Y.cols());
        for (int k = 0; k < d.Y.cols(); ++k) {
            all.push_back(d.Y.col(k));
        }
    }

    Mat3X Y(3, (int)all.size());
    for (int k = 0; k < (int)all.size(); ++k) {
        Y.col(k) = all[k];
    }

    MagData out;
    out.Y = std::move(Y);
    out.source_files = std::move(files);

    if (opt.verbose) {
        std::cout << "Loaded " << out.Y.cols() << " total mag samples from directory " << dirpath
                  << " (files=" << out.source_files.size() << ")\n";
    }
    return out;
}

} // namespace magcal
