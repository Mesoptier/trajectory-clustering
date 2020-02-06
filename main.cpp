#include <iostream>
#include <chrono>
#include <fstream>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"
#include "src/SymmetricMatrix.h"
#include "src/clustering/pam.h"

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

            IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100);
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
    std::cout << "read_curves: " << duration << "ms\n";
}

int main() {
    const auto distance_matrix = read_matrix("data/out/distance_matrix.mtx");
    compute_clusters(distance_matrix, 19); // "characters" has 19 classes

//    const auto curves = read_curves("data/characters/data");
//    const auto distance_matrix = compute_distance_matrix(curves);
//    export_matrix(distance_matrix, "data/out/distance_matrix.mtx");

    return 0;
}