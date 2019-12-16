#include <boost/timer/timer.hpp>
#include <iostream>

#include "io.h"
#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"

int main() {
    boost::timer::auto_cpu_timer timer;

    const Curve curve1 = io::readCurve("data/characters/data/a0001.txt");
    const Curve curve2 = io::readCurve("data/characters/data/a0002.txt");

    std::cout << "reading curves:";
    timer.report();
    timer.start();

    IntegralFrechet alg(curve1, curve2);
    distance_t cost;
    Points matching;
    std::tie(cost, matching) = alg.compute_matching();

    std::cout << "compute matching:";
    timer.report();
    timer.start();

    std::cout << "cost: " << cost << "\n";
    std::cout << "matching size: " << matching.size() << "\n";
//    std::cout << "matching: ";
//    for (auto p : matching) std::cout << "\n| " << p;
//    std::cout << std::endl;

    io::exportPoints("data/out/curve1.csv", curve1.get_points());
    io::exportPoints("data/out/curve2.csv", curve2.get_points());
    io::exportPoints("data/out/matching.csv", matching);

    std::cout << "total:";
    return 0;
}