#include "Curve.h"
#include "FastMarchCDTW.h"

int main() {
    Curve<double> curve1({{0, 0}, {1, 0}});

    Curve<double> curve2({{0.5, 0.5}, {1.5, 0.5}});

    std::cout << FastMarchCDTW::compute(curve1, curve2, 0.05, true) << std::endl;

    return 0;
}