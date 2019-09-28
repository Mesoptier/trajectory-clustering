#include "Curve.h"
#include "FastMarchIntegralFrechet.h"

int main() {
    Curve<double> curve2({{0, 0}, {1, 0}, {1.5, 0.2}});

    Curve<double> curve1({{0, 1}, {2, 2}});

    FastMarchIntegralFrechet solver(curve1, curve2, 0.01, 2, 1);

    std::cout << solver.computeDistance() << std::endl;
    solver.computeMatching();
    solver.save();

    return 0;
}