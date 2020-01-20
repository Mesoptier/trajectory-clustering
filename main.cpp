#include <iostream>
#include <chrono>
#include <fstream>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"

struct MatchingStat {
    std::string name;
    distance_t cost;
    unsigned long long int size;
    long long int duration;
};

int main() {
    //
    // READ CURVES
    //
    auto start_read_curves = std::chrono::high_resolution_clock::now();

    const auto curve1 = io::read_curve("data/characters/data/a0001.txt").coarse();
    const auto curve2 = io::read_curve("data/characters/data/a0002.txt").coarse();

    auto end_read_curves = std::chrono::high_resolution_clock::now();
    std::cout << "read_curves: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_read_curves - start_read_curves).count() << "ms\n";

    //
    // COMPUTE MATCHING
    //

    auto start_compute_matching = std::chrono::high_resolution_clock::now();

    IntegralFrechet alg(curve1, curve2, ParamMetric::L1, 1);
    auto [cost, matching] = alg.compute_matching();

    auto end_compute_matching = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_compute_matching - start_compute_matching).count();

    std::cout << "cost: " << cost << " size: " << matching.size() << " duration: " << duration << "ms\n";

    io::export_points("data/out/curve1.csv", curve1.get_points());
    io::export_points("data/out/curve2.csv", curve2.get_points());
    io::export_points("data/out/matching.csv", matching);

    return 0;
}