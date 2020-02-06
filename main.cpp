#include <iostream>
#include <chrono>
#include <fstream>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"
#include "src/SymmetricMatrix.h"
#include "src/clustering/pam.h"

struct MatchingStat {
    long long int duration;

    // Curves
    std::string curve1_name;
    std::string curve2_name;
    size_t curve1_size;
    size_t curve2_size;

    // Matching
    distance_t cost;
    size_t matching_size;

    // Shortest path search
    size_t nodes_opened;
    size_t nodes_handled;
    size_t nodes_skipped;
};

int main() {
    ///
    // READ DISTANCE MATRIX
    //
    auto start_read_mat = std::chrono::high_resolution_clock::now();

    std::ifstream file("data/out/distance_matrix.mtx");
    auto dissimilarity_matrix = SymmetricMatrix::read(file);

    auto end_read_mat = std::chrono::high_resolution_clock::now();
    std::cout << "read_mat: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_read_mat - start_read_mat).count() << "ms\n";

    //
    // COMPUTE CLUSTERS
    //
    auto start_clustering = std::chrono::high_resolution_clock::now();

    clustering::pam::compute(dissimilarity_matrix.n, 19, dissimilarity_matrix);

    auto end_clustering = std::chrono::high_resolution_clock::now();
    std::cout << "clustering: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_clustering - start_clustering).count() << "ms\n";

    return 0;

//    //
//    // READ CURVES
//    //
//    auto start_read_curves = std::chrono::high_resolution_clock::now();
//
//    const auto curves = io::read_curves("data/characters/data");
//
//    auto end_read_curves = std::chrono::high_resolution_clock::now();
//    std::cout << "read_curves: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_read_curves - start_read_curves).count() << "ms\n";
//
//    //
//    // COMPUTE DISTANCE MATRIX
//    //
//
//    size_t n = curves.size();
//    SymmetricMatrix distance_matrix(n);
//
//    auto start_all = std::chrono::high_resolution_clock::now();
//
//    for (size_t i = 0; i < n; ++i) {
//        const auto curve1 = curves[i].coarse();
//
//        // Don't compute the j == i case, since the cost is always 0 (comparing identical curves)
//        distance_matrix.at(i, i) = 0;
//
//        for (size_t j = i + 1; j < n; ++j) {
//            const auto curve2 = curves[j].coarse();
//
//            IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100);
//            const auto result = alg.compute_matching();
//
//            // TODO: Weigh by max or sum?
//            distance_t averaged_cost = result.cost / std::max(curve1.curve_length(), curve2.curve_length());
//            distance_matrix.at(i, j) = averaged_cost;
//        }
//
//        std::cout << i << '\n';
//    }
//
//    auto end_all = std::chrono::high_resolution_clock::now();
//    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count();
//    std::cout << "total_duration: " << total_duration << "ms\n";
//
//    //
//    // EXPORT DISTANCE MATRIX
//    //
//
//    std::ofstream file("data/out/distance_matrix.mtx");
//    distance_matrix.write(file);
//    file.close();
//
//    return 0;
}