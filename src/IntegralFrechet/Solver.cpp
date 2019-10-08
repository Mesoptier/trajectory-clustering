#include <memory>
#include "Solver.h"

template<class V>
Solver<V>::Solver(const Curve<V>& curve1, const Curve<V>& curve2, double h)
    : curve1(curve1), curve2(curve2), n1(curve1.getNoVertices() - 1), n2(curve2.getNoVertices() - 1) {

    // Create cell grid
    cells.reserve(n1 * n2);

    std::shared_ptr<arma::Col<V>> in1;
    std::shared_ptr<arma::Col<V>> in2;

    for (int i1 = 0; i1 < n1; ++i1) {
        const auto edge1 = curve1.getEdge(i1);
        unsigned int m1 = ceil(edge1.getLength() / h);
        in1 = std::make_shared<arma::Col<V>>(m1);
        in1->fill(INFINITY);
        if (i1 == 0) {
            in1->at(0) = 0;
        }

        for (int i2 = 0; i2 < n2; ++i2) {
            const auto edge2 = curve2.getEdge(i2);
            unsigned int m2 = ceil(edge2.getLength() / h);

            if (i1 == 0) {
                in2 = std::make_shared<arma::Col<V>>(m2);
            } else {
                in2 = cells[(i1 - 1) * n2 + i2].out2;
            }

            const Cell<V> cell(edge1, edge2, m1, m2, in1, in2);
            cells.push_back(cell);

            in1 = cell.out1;

            std::cout << cell << std::endl;
        }
    }
}

template class Solver<double>;