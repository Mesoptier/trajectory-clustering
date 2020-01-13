#include <iostream>
#include <chrono>
#include <fstream>

#include "io.h"
#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"

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

    const auto curves = io::read_curves("data/characters/data");

    auto end_read_curves = std::chrono::high_resolution_clock::now();
    std::cout << "read_curves: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_read_curves - start_read_curves).count() << "ms\n";

    //
    // COMPUTE MATCHING
    //

    const auto curve1 = curves.front().coarse();

    std::vector<MatchingStat> stats;

    for (const auto& curve2 : curves) {
        auto start_compute_matching = std::chrono::high_resolution_clock::now();

        IntegralFrechet alg(curve1, curve2.coarse());
        auto [cost, matching] = alg.compute_matching();

        auto end_compute_matching = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_compute_matching - start_compute_matching).count();

        stats.push_back({
            curve2.name(),
            cost,
            matching.size(),
            duration
        });
    }

    // Save stats to file
    std::ofstream file("data/out/matching_stats.csv");
    file.precision(10);
    for (const auto& stat : stats) {
        file << stat.name << ',' << stat.cost << ',' << stat.size << ',' << stat.duration << '\n';
    }
    file.close();

    return 0;
}