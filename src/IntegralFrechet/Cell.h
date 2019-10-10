#ifndef CODE_CELL_H
#define CODE_CELL_H

#include <utility>
#include <armadillo>
#include <memory>
#include <ostream>
#include "../Edge.h"

template<class V>
class Cell
{
    const Edge<V> edge1;
    const Edge<V> edge2;

    // Number of points "sampled" along the edges (= size of input and output costs)
    int n1;
    int n2;

    // Total costs of getting to an input point (computed by other cells)
    std::shared_ptr<arma::Col<V> const> in1;
    std::shared_ptr<arma::Col<V> const> in2;

public:
    // Total costs of getting to an output point
    std::shared_ptr<arma::Col<V>> out1;
    std::shared_ptr<arma::Col<V>> out2;

    Cell(const Edge<V>& edge1, const Edge<V>& edge2, int n1, int n2, const std::shared_ptr<const arma::Col<V>>& in1,
         const std::shared_ptr<const arma::Col<V>>& in2
    );

    friend std::ostream& operator<<(std::ostream& os, const Cell<V>& cell) {
        return os
            << "----------" << std::endl
            << "Edge 1:\n" << cell.edge1
            << "N: " << cell.n1 << std::endl
            << std::endl
            << "Edge 2:\n" << cell.edge2
            << "N: " << cell.n2 << std::endl
            << std::endl
            << "In 1: " << cell.in1->t()
            << "In 2: " << cell.in2->t()
            << "Out 1:" << cell.out1->t()
            << "Out 2:" << cell.out2->t()
            << std::endl;
    }

private:
    const V& inValue(int i) const;
    V& outValue(int i) const;

    double computeCost(int i, int o) const;
};

#endif //CODE_CELL_H
