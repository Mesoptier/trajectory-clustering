#ifndef CODE_CELL_H
#define CODE_CELL_H

#include <utility>
#include <armadillo>
#include <memory>
#include <ostream>
#include "../Edge.h"
#include "metrics.h"
#include "../geom.h"

class Cell
{
    const Edge edge1;
    const Edge edge2;

    // Number of points "sampled" along the edges (= size of input and output costs)
    int n1;
    int n2;

    // Total costs of getting to an input point (computed by other cells)
    std::shared_ptr<arma::Col<distance_t> const> in1;
    std::shared_ptr<arma::Col<distance_t> const> in2;

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
    std::shared_ptr<arma::Col<distance_t>> out1;
    std::shared_ptr<arma::Col<distance_t>> out2;

    Cell(const Edge& edge1, const Edge& edge2, int n1, int n2, const std::shared_ptr<const arma::Col<distance_t>>& in1,
         const std::shared_ptr<const arma::Col<distance_t>>& in2, Point offset, ImageMetric imageMetric,
         ParamMetric paramMetric
    );

    distance_t getResult() const;
    Points getPath(int i, int o) const;
    Points getMinPath(Point target) const;
    arma::Mat<distance_t> getBoundaryCosts() const;

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
//            Points points;
//            steepestDescent(points, outPoint(i), inPoint(0));
//            writer.writePoints(points);
        }
        writer.closeFunction();
        writer.closeRule();

        // Ellipse axes
        writer.openRule("EllM");
        writer.writeLine(ellipseAxis);
        writer.closeRule();

        writer.openRule("EllH");
        writer.writeLine(ellH);
        writer.closeRule();

        writer.openRule("EllV");
        writer.writeLine(ellV);
        writer.closeRule();

        // Boundary costs



        writer.closeFunction();
    }

    friend std::ostream& operator<<(std::ostream& os, const Cell& cell) {
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
    const distance_t& inValue(int i) const;
    distance_t& outValue(int i) const;

    Point inPoint(int i) const;
    Point outPoint(int i) const;
    int inIndex(const Point& p) const;
    int outIndex(const Point& p) const;

    void steepestDescent(Points& path, Point s, Point t) const;
    distance_t computeCost(int i, int o) const;

    /**
     * Compute the value (cost * distance) of integration along the specified param-edge.
     */
    distance_t integrate(const Point& p1, const Point& p2, ImageMetric imageMetric) const;

    /**
     * Compute the cost (without distance) of linear travel over two edges in image space parametrized by
     * the x- and y- components of the difference vectors at the start/end of the edges.
     */
    distance_t integrate(distance_t dx1, distance_t dy1, distance_t dx2, distance_t dy2, ImageMetric imageMetric) const;
};

#endif //CODE_CELL_H
