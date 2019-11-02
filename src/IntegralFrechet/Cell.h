#ifndef CODE_CELL_H
#define CODE_CELL_H

#include <utility>
#include <armadillo>
#include <memory>
#include <ostream>
#include "../Edge.h"
#include "metrics.h"

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

    arma::Row<V> offset;
    ImageMetric imageMetric;
    ParamMetric paramMetric;

    // Ellipse
    arma::Row<V> midPoint;
    V slope;
    V intercept;

    arma::Col<int> outOrigin;

public:
    // Total costs of getting to an output point
    std::shared_ptr<arma::Col<V>> out1;
    std::shared_ptr<arma::Col<V>> out2;

    Cell(const Edge<V>& edge1, const Edge<V>& edge2, int n1, int n2, const std::shared_ptr<const arma::Col<V>>& in1,
         const std::shared_ptr<const arma::Col<V>>& in2, arma::Row<V> offset, ImageMetric imageMetric,
         ParamMetric paramMetric
    );

    V getResult() const;
    arma::Mat<V> getPath(int i, int o) const;
    arma::Mat<V> getMinPath(arma::Row<V> target) const;

    arma::Row<V> getOffset() {
        return offset;
    }

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

    arma::Row<V> inPoint(int i) const;
    arma::Row<V> outPoint(int i) const;
    int inIndex(arma::Row<V> p) const;
    int outIndex(arma::Row<V> p) const;

    V computeCost(int i, int o) const;

    /**
     * Compute the value of integration along the specified edge.
     *
     * @param p1 Startpoint
     * @param p2 Endpoint
     * @return
     */
    V integrate(arma::Row<V> p1, arma::Row<V> p2) const;
};

#endif //CODE_CELL_H
