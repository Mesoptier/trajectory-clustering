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
        unsigned int m1 = ceil(edge1.length / h) + 1;
        in1 = std::make_shared<arma::Col<V>>(m1);
        in1->fill(INFINITY);
        if (i1 == 0) {
            in1->at(0) = 0;
        }

        for (int i2 = 0; i2 < n2; ++i2) {
            const auto edge2 = curve2.getEdge(i2);
            unsigned int m2 = ceil(edge2.length / h) + 1;

            if (i1 == 0) {
                in2 = std::make_shared<arma::Col<V>>(m2);
                in2->fill(INFINITY);
            } else {
                in2 = cells[(i1 - 1) * n2 + i2].out2;
            }

            arma::Row<V> offset = {curve1.getLength(i1), curve2.getLength(i2)};
            const Cell<V> cell(edge1, edge2, m1, m2, in1, in2, offset);
            cells.push_back(cell);

            in1 = cell.out1;

            std::cout << cell << std::endl;
        }
    }
}

template<class V>
V Solver<V>::getDistance() const  {
    return cells[n1 * n2 - 1].getResult();
}

template<class V>
arma::Mat<V> Solver<V>::getMatching() const {
    arma::Row<V> start = {curve1.getLength(), curve2.getLength()};
    arma::Mat<V> matching(0, 2);

    int i1 = n1 - 1;
    int i2 = n2 - 1;
    bool done = false;

    while (!done) {
        // Do a final loop for the final cell
        done = i1 == 0 && i2 == 0;

        // Get cell
        auto cell = cells[i1 * n2 + i2];

        // Get cell offset
        arma::Row<V> offset = cell.getOffset();

        // Get minimal path to the target point
        auto cellMatching = cell.getMinPath(start - offset);

        // Find index of next cell
        if (cellMatching(0, 0) == 0) {
            i1--;
        }
        if (cellMatching(0, 1) == 0) {
            i2--;
        }

        // Undo translation
        cellMatching.each_row() += offset;

        // Get current endpoint of the matching
        start = cellMatching.row(0);

        // Prepend minimal path to matching
        matching = arma::join_cols(cellMatching, matching);
    }

    return matching;
}

template class Solver<double>;