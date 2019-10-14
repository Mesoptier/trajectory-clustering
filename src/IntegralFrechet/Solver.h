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

    // Grid of cells
    std::vector<Cell<V>> cells;

public:
    Solver(const Curve<V>& curve1, const Curve<V>& curve2, double h);

    V getDistance() const {
        return cells[n1 * n2 - 1].getResult();
    }
};

#endif //CODE_SOLVER_H
