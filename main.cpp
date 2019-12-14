#include <io.h>
#include <Curve.h>
#include <IntegralFrechet/IntegralFrechet.h>
//#include "src/IntegralFrechet/Solver.h"
//#include "FastMarchIntegralFrechet.h"

int main() {
//    const Curve curve1({{0, 0}, {1, 1}, {2, 2}});
//    const Curve curve2({{0, 1}, {1, 1}, {2, 1}});

    const Curve curve1({{0, 0}, {1, 1}, {2, 2}});
    const Curve curve2({{0, 1}, {2, 1}});

    IntegralFrechet alg(curve1, curve2);
    alg.findPath();

//    const Curve curve1 = io::readCurve("data/characters/data/a0001.txt");
//    const Curve curve2 = io::readCurve("data/characters/data/a0002.txt");

//    double h = .1;
//
//    Solver solver(curve1, curve2, h, ImageMetric::L2, ParamMetric::LInfinity_NoShortcuts);
//    auto matching = solver.getMatching();
//
//    std::cout << "Distance: " << solver.getDistance() << std::endl;
//    std::cout << matching << std::endl;

//    io::exportPoints("curve1.csv", curve1.getPoints());
//    io::exportPoints("curve2.csv", curve2.getPoints());
//    io::exportPoints("matching.csv", matching);
//
//    // Solver summary
//    std::ofstream solverFile;
//    solverFile.open("solver.xml");
//
//    ExpressionML::Writer writer(solverFile);
//    writer.open();
//    solver.writeExpressionML(writer);
//    writer.close();
//
//    solverFile.close();


//    double h = 0.01;
//    int imageNorm = 2;
//    int paramNorm = 1;
//
//    FastMarchIntegralFrechet solver(curve1, curve2, h, imageNorm, paramNorm);
//
//    std::cout << "curve1 distance to curve2: " << solver.computeDistance() << std::endl;
//    std::cout << std::endl;
//    solver.computeMatching(0.01);
//    solver.save();
//
//    const int numCenters = 51;
//    std::vector<double> centerRatios;
//    std::vector<Curve<double>> centers;
//
//    for (int i = 0; i < numCenters; ++i) {
//        double ratio = (double) i / (numCenters - 1);
//        solver.computeCenter(ratio);
//        centers.push_back(Curve<double>(solver.getCenter()));
//        centerRatios.push_back(ratio);
//    }
//
//    arma::mat centerDistances(centers.size(), 2);
//
//    double minSum = INFINITY;
//    double minSumRatio = 0;
//    Curve<double> minSumCenter({});
//
//    double minDiff = INFINITY;
//    double minDiffRatio = 0;
//    Curve<double> minDiffCenter({});
//
//    double minSumSquares = INFINITY;
//    double minSumSquaresRatio = 0;
//    Curve<double> minSumSquaresCenter({});
//
//    for (int i = 0; i < centers.size(); ++i) {
//        const auto& center = centers[i];
//
//        solver = FastMarchIntegralFrechet(curve1, center, h, imageNorm, paramNorm);
//        centerDistances(i, 0) = solver.computeDistance();
//
//        solver = FastMarchIntegralFrechet(curve2, center, h, imageNorm, paramNorm);
//        centerDistances(i, 1) = solver.computeDistance();
//
//        double sum = centerDistances(i, 0) + centerDistances(i, 1);
//        if (sum < minSum) {
//            minSum = sum;
//            minSumCenter = center;
//            minSumRatio = centerRatios[i];
//        }
//
//        double diff = std::abs(centerDistances(i, 0) - centerDistances(i, 1));
//        if (diff < minDiff) {
//            minDiff = diff;
//            minDiffCenter = center;
//            minDiffRatio = centerRatios[i];
//        }
//
//        double sumSquares = pow(centerDistances(i, 0), 2) + pow(centerDistances(i, 1), 2);
//        if (sumSquares < minSumSquares) {
//            minSumSquares = sumSquares;
//            minSumSquaresCenter = center;
//            minSumSquaresRatio = centerRatios[i];
//        }
//    }
//
//    arma::rowvec({minSumRatio, minDiffRatio, minSumSquaresRatio}).save("minCenterValues.csv",arma::csv_ascii);
//    minSumCenter.getVertices().save("minSumCenter.csv", arma::csv_ascii);
//    minDiffCenter.getVertices().save("minDiffCenter.csv", arma::csv_ascii);
//    minSumSquaresCenter.getVertices().save("minSumSquaresCenter.csv", arma::csv_ascii);
//    centerDistances.save("centerDistances.csv", arma::csv_ascii);

    return 0;
}