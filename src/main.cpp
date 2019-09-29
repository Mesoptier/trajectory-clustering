#include "Curve.h"
#include "FastMarchIntegralFrechet.h"

int main() {
    const Curve<double> curve1({{0, 1}, {2, 3}});
    const Curve<double> curve2({{0, 0}, {1, 0}, {1.5, 0.2}});

    double h = 0.01;
    int imageNorm = 2;
    int paramNorm = 1;

    FastMarchIntegralFrechet solver(curve1, curve2, h, imageNorm, paramNorm);

    std::cout << "curve1 distance to curve2: " << solver.computeDistance() << std::endl;
    std::cout << std::endl;
    solver.computeMatching(0.1);
    solver.save();

    const int numCenters = 51;
    std::vector<double> centerRatios;
    std::vector<Curve<double>> centers;

    for (int i = 0; i < numCenters; ++i) {
        double ratio = (double) i / (numCenters - 1);
        solver.computeCenter(ratio);
        centers.push_back(Curve<double>(solver.getCenter()));
        centerRatios.push_back(ratio);
    }

    arma::mat centerDistances(centers.size(), 2);

    double minSum = INFINITY;
    double minSumRatio = 0;
    Curve<double> minSumCenter({});

    double minDiff = INFINITY;
    double minDiffRatio = 0;
    Curve<double> minDiffCenter({});

    for (int i = 0; i < centers.size(); ++i) {
        const auto& center = centers[i];

        solver = FastMarchIntegralFrechet(curve1, center, h, imageNorm, paramNorm);
        centerDistances(i, 0) = solver.computeDistance();

        solver = FastMarchIntegralFrechet(curve2, center, h, imageNorm, paramNorm);
        centerDistances(i, 1) = solver.computeDistance();

        double sum = centerDistances(i, 0) + centerDistances(i, 1);
        if (sum < minSum) {
            minSum = sum;
            minSumCenter = center;
            minSumRatio = centerRatios[i];
        }

        double diff = std::abs(centerDistances(i, 0) - centerDistances(i, 1));
        if (diff < minDiff) {
            minDiff = diff;
            minDiffCenter = center;
            minDiffRatio = centerRatios[i];
        }
    }

    arma::rowvec({minSumRatio, minDiffRatio}).save("minCenterValues.csv",arma::csv_ascii);
    minSumCenter.getVertices().save("minSumCenter.csv", arma::csv_ascii);
    minDiffCenter.getVertices().save("minDiffCenter.csv", arma::csv_ascii);
    centerDistances.save("centerDistances.csv", arma::csv_ascii);

    return 0;
}