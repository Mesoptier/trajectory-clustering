#include <iostream>
#include <chrono>
#include <fstream>

#include "src/io.h"
#include "src/Curve.h"
#include "src/IntegralFrechet/IntegralFrechet.h"

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

    const auto& curve1 = curves.front();

    std::vector<MatchingStat> stats;

    auto start_all = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < curves.size(); ++i) {
        const auto curve2 = curves[i];

        auto start_compute_matching = std::chrono::high_resolution_clock::now();

        IntegralFrechet alg(curve1, curve2, ParamMetric::LInfinity_NoShortcuts, 100);
        const auto result = alg.compute_matching();

        auto end_compute_matching = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_compute_matching - start_compute_matching).count();

        stats.push_back({
            duration,

            curve1.name(),
            curve2.name(),
            curve1.size(),
            curve2.size(),

            result.cost,
            result.matching.size(),

            result.search_stat.nodes_opened,
            result.search_stat.nodes_handled,
            result.search_stat.nodes_skipped,
        });
    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count();
    std::cout << "total_duration: " << total_duration << "ms\n";

    //
    // EXPORT MATCHING STATS
    //

    std::ofstream file("data/out/matching_stats.csv");

    file.precision(10);
    for (const auto& stat : stats) {
        file << stat.curve1_size << ',' << stat.curve2_size << ',' << stat.nodes_opened << '\n';
    }

    file.close();

    return 0;
}