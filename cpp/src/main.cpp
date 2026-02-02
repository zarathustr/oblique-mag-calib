#include "magcal/magcal.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace magcal;

struct Args {
    std::string input;
    std::string outdir = "out";
    int downsample = 1;
    double lambda = 0.0;
    int maxiter = 200;
    bool only_init = false;
    int nu = 60;
    int nv = 30;
    bool verbose = true;
};

static void print_usage(const char* prog) {
    std::cout
        << "Usage: " << prog << " --input <csv_txt_or_dir> [options]\n\n"
        << "Options:\n"
        << "  --outdir <path>      Output directory (default: out)\n"
        << "  --downsample <k>     Keep every k-th sample (default: 1)\n"
        << "  --lambda <val>       Smoothness penalty on M (default: 0)\n"
        << "  --maxiter <k>        Max iterations for manifold solver (default: 200)\n"
        << "  --only_init          Skip manifold refinement (ellipsoid init only)\n"
        << "  --nu <k>             Ellipsoid mesh U resolution (default: 60)\n"
        << "  --nv <k>             Ellipsoid mesh V resolution (default: 30)\n"
        << "  --quiet              Less console output\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](const std::string& name) {
            if (i + 1 >= argc) throw std::runtime_error("Missing value after " + name);
            return std::string(argv[++i]);
        };
        if (s == "--input") a.input = need(s);
        else if (s == "--outdir") a.outdir = need(s);
        else if (s == "--downsample") a.downsample = std::stoi(need(s));
        else if (s == "--lambda") a.lambda = std::stod(need(s));
        else if (s == "--maxiter") a.maxiter = std::stoi(need(s));
        else if (s == "--nu") a.nu = std::stoi(need(s));
        else if (s == "--nv") a.nv = std::stoi(need(s));
        else if (s == "--only_init") a.only_init = true;
        else if (s == "--quiet") a.verbose = false;
        else if (s == "--help" || s == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown option: " + s);
        }
    }
    if (a.input.empty()) {
        throw std::runtime_error("--input is required");
    }
    if (a.downsample < 1) a.downsample = 1;
    if (a.maxiter < 1) a.maxiter = 1;
    return a;
}

int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);

        std::filesystem::create_directories(args.outdir);

        CsvLoadOptions lopt;
        lopt.downsample = args.downsample;
        lopt.verbose = args.verbose;

        Txt10LoadOptions topt;
        topt.downsample = args.downsample;
        topt.verbose = args.verbose;

        MagData data;
        if (std::filesystem::is_directory(args.input)) {
            // Directory may contain a mix of .csv and .txt logs.
            data = load_magnetometer_directory_mixed(
                args.input,
                [&](const std::string& f) { return load_magnetometer_csv(f, lopt); },
                topt);
        } else {
            // Auto-detect by file extension.
            const auto ext = to_lower(std::filesystem::path(args.input).extension().string());
            if (ext == ".txt") {
                data = load_nautilus_txt10(args.input, topt);
            } else {
                data = load_magnetometer_csv(args.input, lopt);
            }
        }

        const Mat3X& Y = data.Y;
        std::cout << "N=" << Y.cols() << " samples\n";

        // Always compute Wu init and export initial ellipsoid
        WuInitResult init = wu_initial_estimate(Y, args.verbose);
        write_matrix_txt(init.T, args.outdir + "/T_init.txt");
        write_matrix_txt(init.A, args.outdir + "/A_init.txt");
        write_vector_txt(init.h, args.outdir + "/h_init.txt");
        write_obj_ellipsoid(init.A, init.h, args.outdir + "/ellipsoid_init.obj", args.nu, args.nv);

        Mat3 T_hat = init.T;
        Vec3 h_hat = init.h;

        if (!args.only_init) {
            ManifoldOptions opt;
            opt.max_iters = args.maxiter;
            opt.lambda_smooth = args.lambda;
            opt.verbose = args.verbose;

            ManifoldResult res = solve_product_manifold(Y, opt);
            T_hat = res.T;
            h_hat = res.h;

            std::cout << "Converged: " << (res.converged ? "yes" : "no")
                      << "  iters=" << res.iters
                      << "  final cost=" << res.cost << "\n";
        }

        // Calibrated unit vectors
        Mat3X Mcal = T_hat.inverse() * (Y - h_hat * Eigen::RowVectorXd::Ones(Y.cols()));

        // Save
        write_matrix_txt(T_hat, args.outdir + "/T_hat.txt");
        write_vector_txt(h_hat, args.outdir + "/h_hat.txt");
        write_calibrated_csv(Mcal, args.outdir + "/calibrated.csv");

        // Export refined ellipsoid mesh
        Mat3 A_hat = ellipsoid_A_from_T(T_hat);
        write_matrix_txt(A_hat, args.outdir + "/A_hat.txt");
        write_obj_ellipsoid(A_hat, h_hat, args.outdir + "/ellipsoid_hat.obj", args.nu, args.nv);

        // Basic diagnostics
        Eigen::RowVectorXd mag_raw = Y.colwise().norm();
        Eigen::RowVectorXd mag_cal = Mcal.colwise().norm();

        double mean_abs_mag_err = (mag_cal.array() - 1.0).abs().mean();
        std::cout << "Mean | ||m_cal|| - 1 | = " << mean_abs_mag_err << "\n";
        std::cout << "Raw magnitude mean/std: " << mag_raw.mean() << " / "
                  << std::sqrt((mag_raw.array() - mag_raw.mean()).square().mean()) << "\n";

        std::cout << "Wrote results to: " << args.outdir << "\n";
        std::cout << "OBJ meshes: ellipsoid_init.obj, ellipsoid_hat.obj\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }
}
