#ifndef CODE_SOLVER_H
#define CODE_SOLVER_H

#include <vector>
#include "Cell.h"
#include "../Curve.h"

template<class V>
class Solver
{
    Curve<double> curve1;
    Curve<double> curve2;

    // Curve complexity
    unsigned int n1;
    unsigned int n2;

    ImageMetric imageMetric;
    ParamMetric paramMetric;

    // Grid of cells
    std::vector<Cell<V>> cells;

public:
    Solver(const Curve<V>& curve1, const Curve<V>& curve2, double h, ImageMetric imageMetric, ParamMetric paramMetric);

    V getDistance() const;
    arma::Mat<V> getMatching() const;
    arma::Mat<V> getBoundaryCosts() const;
};

#endif //CODE_SOLVER_H
