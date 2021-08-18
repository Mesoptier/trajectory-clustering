#include <iostream>
#include <chrono>
#include <fstream>
#include <random>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"
#include "src/SymmetricMatrix.h"
#include "src/clustering/pam.h"
#include "src/IntegralFrechet/MatchingBand.h"
#include "src/simplification/imaiiri.h"
#include "src/DTW/dtw.h"
#include "src/cdtw/cdtw.h"
#include "src/cdtw/1d-l1-l1.h"
#include "src/cdtw/2d-l1-l1.h"
#include "src/cdtw/1d-l2squared-l1.h"
#include "src/cdtw/2d-l2squared-l1.h"

namespace {
//
// I/O helpers
//

std::vector<Curve> read_curves(const std::string& filename) {
    auto start_time = std::chrono::high_resolution_clock::now();

    const auto curves = io::read_curves(filename);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "read_curves: " << duration << "ms\n";

    return curves;
}

SymmetricMatrix read_matrix(const std::string& filename) {
    auto start_time = std::chrono::high_resolution_clock::now();

    std::ifstream file(filename);
    const auto matrix = SymmetricMatrix::read(file);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "read_matrix: " << duration << "ms\n";

    return matrix;
}

void export_matrix(const SymmetricMatrix& matrix, const std::string& filename) {
    std::ofstream file(filename);
    matrix.write(file);
    file.close();
}


//
// Algorithms
//

SymmetricMatrix compute_distance_matrix(const std::vector<Curve>& curves) {
    size_t n = curves.size();
    SymmetricMatrix distance_matrix(n);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < n; ++i) {
        const auto curve1 = curves[i];

        // Don't compute the j == i case, since the cost is always 0 (comparing identical curves)
        distance_matrix.at(i, i) = 0;

        for (size_t j = i + 1; j < n; ++j) {
            const auto curve2 = curves[j];

            IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100, nullptr);
            const auto result = alg.compute_matching();

            // TODO: Weigh by max or sum?
            // TODO: Move to IntegralFrechet class
            distance_t averaged_cost = result.cost / std::max(curve1.curve_length(), curve2.curve_length());
            distance_matrix.at(i, j) = averaged_cost;
        }

        std::cout << i << '\n';
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_distance_matrix: " << duration << "ms\n";

    return distance_matrix;
}

void compute_clusters(const SymmetricMatrix& distance_matrix, size_t k) {
    auto start_time = std::chrono::high_resolution_clock::now();

    clustering::pam::compute(distance_matrix.n, k, distance_matrix);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_clusters: " << duration << "ms\n";
}

IntegralFrechet::MatchingResult compute_matching(const Curve& curve1, const Curve& curve2, distance_t resolution, const MatchingBand* const band) {
    auto start_time = std::chrono::high_resolution_clock::now();

    IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, resolution, band);
    const auto result = alg.compute_matching();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_matching: " << duration << "ms\n";

    return result;
}

MatchingBand compute_band(const Curve& curve1, const Curve& curve2, const Points& matching, distance_t radius) {
    auto start_time = std::chrono::high_resolution_clock::now();

    const MatchingBand band(curve1, curve2, matching, radius);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "compute_band: " << duration << "ms\n";

    return band;
}

//
// Experiments
//

void experiment_with_or_without_bands() {
    struct Stat {
        distance_t band_cost;
        long long int band_time;
        distance_t no_band_cost;
        long long int no_band_time;
    };
    std::vector<Stat> stats;

    bool maintain_lengths = true;
    distance_t simple_resolution = 10;
    distance_t resolution = 1;
    distance_t band_radius = 1;

    const auto curves = read_curves("data/characters/data");

    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_int_distribution<size_t> distribution(0, curves.size() - 1);

    size_t k_max = 1000;
    for (size_t k = 0; k < k_max; ++k) {
        auto i = distribution(generator);
        auto j = distribution(generator);

        const auto curve1 = curves.at(i);
        const auto curve2 = curves.at(j);

        // With band
        auto start_band = std::chrono::high_resolution_clock::now();
        const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), simple_resolution,nullptr);
        const auto band = compute_band(curve1, curve2, result_alt.matching, band_radius);
        const auto result_band = compute_matching(curve1, curve2, resolution, &band);
        auto end_band = std::chrono::high_resolution_clock::now();
        auto duration_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_band - start_band).count();

        // Without band
        auto start_no_band = std::chrono::high_resolution_clock::now();
        const auto result_no_band = compute_matching(curve1, curve2, resolution, nullptr);
        auto end_no_band = std::chrono::high_resolution_clock::now();
        auto duration_no_band = std::chrono::duration_cast<std::chrono::milliseconds>(end_no_band - start_no_band).count();

        stats.push_back({
            result_band.cost,
            duration_band,
            result_no_band.cost,
            duration_no_band,
        });

        std::cout << (k + 1) << " / " << k_max << '\n';
    }

    std::ofstream file("data/out/band_stats.tsv");
    file << "band_cost\tband_time\tno_band_cost\tno_band_time\n";
    for (const auto& stat : stats) {
        file << stat.band_cost << '\t' << stat.band_time << '\t' << stat.no_band_cost << '\t' << stat.no_band_time << '\n';
    }
    file.close();
}

void experiment_visualize_band() {
    bool maintain_lengths = true;

    const auto curve1 = io::read_curve("data/characters/data/a0001.txt");
    const auto curve2 = io::read_curve("data/characters/data/a0002.txt");

    const auto result_no_band = compute_matching(curve1, curve2, 1, nullptr);
    io::export_points("data/out/debug_points_old.csv", result_no_band.search_stat.nodes_as_points);
    io::export_points("data/out/matching3.csv", result_no_band.matching);
    std::cout << "matching (no band) cost: " << result_no_band.cost << '\n';

    auto start_time = std::chrono::high_resolution_clock::now();

    const auto result_alt = compute_matching(curve1.simplify(maintain_lengths), curve2.simplify(maintain_lengths), 10,nullptr);
    const auto band = compute_band(curve1, curve2, result_alt.matching,     2);
    const auto result = compute_matching(curve1, curve2, 1, &band);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "main: " << duration << "ms\n";

    io::export_points("data/out/curve1.csv", curve1.get_points());
    io::export_points("data/out/curve2.csv", curve2.get_points());
    io::export_points("data/out/matching1.csv", result.matching);
    io::export_points("data/out/matching2.csv", result_alt.matching);
    io::export_points("data/out/debug_points.csv", result.search_stat.nodes_as_points);

    std::cout << "matching cost: " << result.cost << '\n';
    std::cout << "matching (alt) cost: " << result_alt.cost << '\n';
}

void experiment_compare_heuristic_vs_extact_cdtw() {
//    const std::string data_dir = "D:/TUe/MSc Project/code/data";
    const std::string data_dir = "data";

    auto curve1 = io::read_curve(data_dir + "/characters/data/a0001.txt");
    auto curve2 = io::read_curve(data_dir + "/characters/data/a0002.txt");

//    const auto curve1 = io::read_curve("data/characters/data/a0003.txt").simplify(false);
//    const auto curve2 = io::read_curve("data/characters/data/a0002.txt").simplify(false);

//    const Curve curve2("segment", {{0, 0}, {10, 10}});

//    io::export_points("data/out/curve1.csv", curve1.get_points());
//    io::export_points("data/out/curve2.csv", curve2.get_points());


    // Heuristic
    {
        auto start_time = std::chrono::high_resolution_clock::now();

        IntegralFrechet heuristic_alg(curve1, curve2, ParamMetric::L1, .1);
        const auto heuristic_res = heuristic_alg.compute_matching();

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "heuristic cost: " << heuristic_res.cost << '\n';
        std::cout << "heuristic time: " << duration << "ms\n";
        std::cout << "heuristic nodes_handled: " << heuristic_res.search_stat.nodes_handled << "\n";
        std::cout << "heuristic nodes_skipped: " << heuristic_res.search_stat.nodes_skipped << "\n";
        std::cout << "heuristic nodes_opened: " << heuristic_res.search_stat.nodes_opened << "\n";
    }

    // Exact
    {
        auto start_time = std::chrono::high_resolution_clock::now();

        CDTW<2, Norm::L2Squared, Norm::L1> exact_alg(curve1, curve2);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//        std::cout << "exact cost: " << exact_alg.cost() << '\n';
        std::cout << "exact time: " << duration << "ms\n";

//        exact_alg.output_visualization_data();
//        exact_alg.print_complexity();
    }

}
}

int main() {
    // TODO: Compare Heuristic CDTW vs Exact CDTW (timing and result)
    // TODO: Upgrade to 2D + L2^2
    // TODO: How to Dijkstra in Exact CDTW?
    auto curve1 = io::read_curve("data/characters/data/a0015.txt");
    auto curve2 = io::read_curve("data/characters/data/b0131.txt");

    Point p1(0, 0);
    Point p2(1, 1);
    Point p3(2, 0);
    Point p4(3, 1);
    Point p5(4, 0);
    Point p6(5, 1);
    Point p7(6, 0);
    Point p8(7, 1);

    Point p9(0, 0.7);
    Point p10(5, 0.5);
    // Point p11(2, 1);
    // Point p12(3, 0);
    // Point p13(4, 1);
    // Point p14(5, 0);
    // Point p15(6, 1);
    // Point p16(7, 0);


    Curve c1("", {p1, p2, p3, p4, p5, p6, p7, p8});
    Curve c2("", {p9, p10/*, p11, p12, p13, p14, p15, p16*/});
    
    using _CDTW = CDTW<2, Norm::L1, Norm::L1>;
    auto cdtw = _CDTW(curve1, curve2);
    // auto cdtw = _CDTW(c1, c2);

    // experiment_compare_heuristic_vs_extact_cdtw();

    return 0;
}
