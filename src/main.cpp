#include "Curve.h"
#include "FastMarchCDTW.h"

int main() {
    Curve<double> curve2({{0, 0}, {1, 0}, {1.5, 0.2}});

    Curve<double> curve1({{0, 1}, {2, 2}});

    FastMarchCDTW solver;

    std::cout << solver.compute(curve1, curve2, 0.01, 2, 1, true) << std::endl;

    return 0;
}