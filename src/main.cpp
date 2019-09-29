#include "Curve.h"
#include "FastMarchIntegralFrechet.h"

int main() {
    Curve<double> curve1({{0, 1}, {2, 3}});
    Curve<double> curve2({{0, 0}, {1, 0}, {1.5, 0.2}});

    double h = 0.01;
    int imageNorm = 2;
    int paramNorm = 1;

    FastMarchIntegralFrechet solver(curve1, curve2, h, imageNorm, paramNorm);

    std::cout << "curve1 distance to curve2: " << solver.computeDistance() << std::endl;
    solver.computeMatching();
    solver.computeCenter();
    solver.save();

    Curve<double> center(solver.getCenter());

    solver = FastMarchIntegralFrechet(curve1, center, h, imageNorm, paramNorm);
    std::cout << "curve1 distance to center: " << solver.computeDistance() << std::endl;

    solver = FastMarchIntegralFrechet(curve2, center, h, imageNorm, paramNorm);
    std::cout << "curve2 distance to center: " << solver.computeDistance() << std::endl;

    return 0;
}