#include <iostream>

#include "io.h"
#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"

int main() {
    const Curve curve1 = io::readCurve("data/characters/data/a0001.txt").coarse();
    const Curve curve2 = io::readCurve("data/characters/data/a0002.txt").coarse();

    IntegralFrechet alg(curve1, curve2);
    distance_t cost;
    Points matching;
    std::tie(cost, matching) = alg.compute_matching();

    std::cout << "cost: " << cost << "\n";
    std::cout << "matching size: " << matching.size() << "\n";
//    std::cout << "matching: ";
//    for (auto p : matching) std::cout << "\n| " << p;
//    std::cout << std::endl;

    io::exportPoints("data/out/curve1.csv", curve1.get_points());
    io::exportPoints("data/out/curve2.csv", curve2.get_points());
    io::exportPoints("data/out/matching.csv", matching);

    return 0;
}