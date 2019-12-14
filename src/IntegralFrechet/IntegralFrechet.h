#pragma once

#include <Curve.h>
#include <geom.h>
#include <ostream>
#include <utility>
#include "metrics.h"

struct Cell {
    // Edges (edge1 = {p1, p1 + 1}, edge2 = {p2, p2 + 1})
//    PointID p1;
//    PointID p2;

//    // Edge lengths
//    distance_t len1;
//    distance_t len2;

    size_t n1;
    size_t n2;

    Line ellM; // Monotone axis
    Line ellH; // Line through points where ellipse tangent is horizontal
    Line ellV; // Line through points where ellipse tangent is vertical

    Cell(const Point& s1, const Point& s2, const Point& t1, const Point& t2) {
        // TODO: Make resolution a parameter
        distance_t resolution = 0.2;
        n1 = ceil(s1.dist(t1) / resolution) + 1;
        n2 = ceil(s2.dist(t2) / resolution) + 1;
    }
};

class IntegralFrechet
{
private:
    const Curve curve1;
    const Curve curve2;

    std::vector<Cell> cells;

    const Cell& get_cell(PointID p1, PointID p2) const;
    Cell get_cell(CPoint p1, CPoint p2) const {
        return get_cell(p1.getPoint(), p2.getPoint());
    };
    Cell get_cell(CPosition pos) const {
        return get_cell(pos[0], pos[1]);
    };

    /**
     * Compute the cost of the optimal path from s through t through a single cell.
     *
     * @tparam imageMetric Which metric to use in image space.
     * @tparam paramMetric Which metric to use in parameter space.
     * @param s Source point.
     * @param t Target point. Should be monotonically greater than s.
     * @return
     */
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t cost(const CPosition& s, const CPosition& t) const;

    // Optimal path from s to t
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    CPositions compute_path(const CPosition& s, const CPosition& t) const;

    // Steepest descent from s in the direction of t
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    CPositions steepest_descent(const CPosition& s, const CPosition& t) const;

    // Compute cost over path by integration
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t integrate(const CPositions& path) const;

    // Cost from distance between curves
    template<ImageMetric imageMetric>
    distance_t integrate_cost(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const;

    // Cost from distance travelled over curve
    template<ParamMetric paramMetric>
    distance_t integrate_dist(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const;

public:

    IntegralFrechet(const Curve& curve1, const Curve& curve2);

    CPositions compute_matching();

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(const Node& node, std::vector<Node>& neighbors) const;
    cost_t cost(const Node& s, const Node& t) const;
    cost_t heuristic_cost(const Node& s, const Node& goal) const;
};