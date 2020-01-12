#include <iostream>
#include <chrono>

#include "io.h"
#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"

int main() {
//    const auto curves = io::read_curves("data/characters/data");
//    size_t count = 0;
//
//    auto start = std::chrono::high_resolution_clock::now();
//
//    for (size_t i = 0; i < curves.size(); ++i) {
//        auto curve1 = curves[i];
//
//        for (size_t j = i + 1; j < curves.size(); ++j) {
//            auto curve2 = curves[j];
//
////            std::cout << "complexity: " << curve1.size() << ' ' << curve2.size() << '\n';
//
//            IntegralFrechet alg(curve1, curve2);
//
//            const auto [cost, matching] = alg.compute_matching();
//
//            if (count % 100 == 0) {
//                auto time = std::chrono::high_resolution_clock::now();;
//                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time - start).count();
//                std::cout << count << ' ' << duration << '\n';
//            }
//
//            ++count;
//        }
//    }
//
//    std::cout << "count: " << count << std::endl;


    //
    // READ CURVES
    //
    auto start_read_curves = std::chrono::high_resolution_clock::now();

    const Curve curve1 = io::read_curve("data/characters/data/a0001.txt");
    const Curve curve2 = io::read_curve("data/characters/data/a0002.txt");

    auto end_read_curves = std::chrono::high_resolution_clock::now();
    std::cout << "read_curves: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_read_curves - start_read_curves).count() << "ms\n";

    //
    // COMPUTE MATCHING
    //
    auto start_compute_matching = std::chrono::high_resolution_clock::now();

    IntegralFrechet alg(curve1, curve2);
    distance_t cost;
    Points matching;
    std::tie(cost, matching) = alg.compute_matching();

    auto end_compute_matching = std::chrono::high_resolution_clock::now();
    std::cout << "compute_matching: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_compute_matching - start_compute_matching).count() << "ms\n";

    std::cout << "cost: " << cost << "\n";
    std::cout << "matching size: " << matching.size() << "\n";
//    std::cout << "matching: ";
//    for (auto p : matching) std::cout << "\n| " << p;
//    std::cout << std::endl;

    io::export_points("data/out/curve1.csv", curve1.get_points());
    io::export_points("data/out/curve2.csv", curve2.get_points());
    io::export_points("data/out/matching.csv", matching);

    return 0;
}