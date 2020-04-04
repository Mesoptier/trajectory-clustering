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
#include "src/cdtw/1d-l1.h"

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
}

int main() {
//    const auto curves = read_curves("data/characters/data");
//    // const auto dm = compute_distance_matrix(curves);
//    // export_matrix(dm, "data/out/distance_matrix.mtx");
//
//    // const auto distance_matrix = read_matrix("data/out/distance_matrix.mtx");
//    // compute_clusters(distance_matrix, 19); // "characters" has 19 classes
//{
//    auto start_time = std::chrono::high_resolution_clock::now();
//    auto ret = simplification::imai_iri::simplify(curves[0], 5);
//    auto end_time = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//    std::cout << "Time: " << duration << "ms\n";
//    std::cout << "Cost: " << ret.first << "\n";
//    std::cout << "Curve:";
//    for (const auto& p: ret.second.get_points())
//        std::cout << " " << p;
//    std::cout << "\n";
//    std::cout << "Length: " << ret.second.size() << std::endl;
//}
//{
//    auto start_time = std::chrono::high_resolution_clock::now();
//    DTW distance(curves[0], curves[1]);
//    auto end_time = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//    std::cout << "Time: " << duration << "ms\n";
//    std::cout << "Cost: " << distance.cost() << "\n";
//    std::cout << "Matching:";
//    for (const auto& p: distance.matching())
//        std::cout << " (" << p.first << ", " << p.second << ")";
//    std::cout << std::endl;
//}

    // To simplify the expressions we assume that (0, 0) in the infinite parameter space corresponds to the "center" of
    // the ellipses (though in 1D these ellipses are degenerate). The cell is then translated in this space to have the
    // correct height function.

    Curve curve1("curve1", {{0, 0}, {6, 0}});
    Curve curve2("curve2", {{0, 0}, {8, 0}, {6, 0}});
    CDTW<1, Norm::L1, Norm::L1> cdtw(curve1, curve2);
//
//    // Bottom-left corner of cell
//    double sx = -5;
//    double sy = 0;
//    // Top-right corner of cell
//    double tx = 5;
//    double ty = 10;
//
//    std::cout << "RIGHT = {\n";
//    std::cout << "\t" << bottom_to_right_1(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_right_2(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_right_3(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_right_4(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_right_5(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_right_6(sx, sy, tx, ty) << ",\n";
//    std::cout << "} -> " << bottom_to_right(sx, sy, tx, ty) << "\n";
//
//    std::cout << "TOP = {\n";
//    std::cout << "\t" << bottom_to_top_1(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_top_2(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_top_3(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_top_4(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_top_5(sx, sy, tx, ty) << ",\n";
//    std::cout << "\t" << bottom_to_top_6(sx, sy, tx, ty) << ",\n";
//    std::cout << "} -> " << bottom_to_top(sx, sy, tx, ty) << "\n";

    return 0;
}
