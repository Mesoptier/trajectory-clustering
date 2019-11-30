#ifndef CODE_CELL_H
#define CODE_CELL_H

#include <utility>
#include <armadillo>
#include <memory>
#include <ostream>
#include "../Edge.h"
#include "metrics.h"
#include "../geom.h"

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

    Point offset;
    ImageMetric imageMetric;
    ParamMetric paramMetric;

    // Ellipse
    Point midPoint;

    Line ellipseAxis; // Monotone ellipse axis
    Line ellH; // Line through points where tangent to ellipse is horizontal
    Line ellV; // Line through points where tangent to ellipse is vertical

    arma::Col<int> outOrigin;

public:
    // Total costs of getting to an output point
    std::shared_ptr<arma::Col<V>> out1;
    std::shared_ptr<arma::Col<V>> out2;

    Cell(const Edge<V>& edge1, const Edge<V>& edge2, int n1, int n2, const std::shared_ptr<const arma::Col<V>>& in1,
         const std::shared_ptr<const arma::Col<V>>& in2, Point offset, ImageMetric imageMetric,
         ParamMetric paramMetric
    );

    V getResult() const;
    arma::Mat<V> getPath(int i, int o) const;
    Points getMinPath(Point target) const;
    arma::Mat<V> getBoundaryCosts() const;

    Point getOffset() const {
        return offset;
    }

    void writeExpressionML(ExpressionML::Writer& writer) const {
        writer.openFunction("Association");

        // Offset
        writer.openRule("Offset");
        writer.writePoint(getOffset());
        writer.closeRule();

        // Ellipse midpoint
        writer.openRule("EllipseMidpoint");
        writer.writePoint(midPoint);
        writer.closeRule();

        // Edges
        writer.openRule("Edges");
        writer.openFunction("List");
        edge1.writeExpressionML(writer);
        edge2.writeExpressionML(writer);
        writer.closeFunction();
        writer.closeRule();

        // Optimal paths (to each sampled point on out-boundary)
        writer.openRule("Paths");
        writer.openFunction("List");
        for (int i = 0; i < n1 + n2; ++i) {
            writer.writePoints(getMinPath(outPoint(i)));
        }
        writer.closeFunction();
        writer.closeRule();

        // Ellipse axes
        writer.openRule("MonotoneEllipseAxis");
        writer.writeLine(ellipseAxis);
        writer.closeRule();

        // Boundary costs



        writer.closeFunction();
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
     * Compute the value (cost * distance) of integration along the specified param-edge.
     */
    V integrate(arma::Row<V> p1, arma::Row<V> p2, ImageMetric imageMetric) const;

    /**
     * Compute the cost (without distance) of linear travel over two edges in image space parametrized by
     * the x- and y- components of the difference vectors at the start/end of the edges.
     */
    V integrate(V dx1, V dy1, V dx2, V dy2, ImageMetric imageMetric) const;
};

#endif //CODE_CELL_H
