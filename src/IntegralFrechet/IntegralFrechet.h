#pragma once

#include <Curve.h>
#include <geom.h>
#include <ostream>
#include <utility>
#include "metrics.h"

using CellCoordinate = std::array<PointID, 2>;

struct Cell {
    // Edges (edge1 = {p1, p1 + 1}, edge2 = {p2, p2 + 1})
//    PointID p1;
//    PointID p2;

//    // Edge lengths
//    distance_t len1;
//    distance_t len2;

    size_t n1;
    size_t n2;

    Point center;
    Line ell_m; // Monotone axis
//    Line ellH; // Line through points where ellipse tangent is horizontal
//    Line ellV; // Line through points where ellipse tangent is vertical

    Cell(const Point& s1, const Point& s2, const Point& t1, const Point& t2) {
        // TODO: Make resolution a parameter
        distance_t resolution = 0.2;
        n1 = ceil(s1.dist(t1) / resolution) + 1;
        n2 = ceil(s2.dist(t2) / resolution) + 1;

        Edge e1(s1, t1);
        Edge e2(s2, t2);

        if (isParallel(e1.line, e2.line)) {
            const auto image_point = e1.line.closest(s2);
            center = {e1.param(image_point), 0};
        } else {
            const auto image_point = intersect(e1.line, e2.line);
            center = {e1.param(image_point), e2.param(image_point)};
        }

        ell_m = {center, {1, 1}};
    }
};

class IntegralFrechet
{
private:
    const Curve curve1;
    const Curve curve2;

    std::vector<Cell> cells;

    const Cell& get_cell(CellCoordinate cc) const;

    /**
     * Convert a global CPosition to a local cell point.
     */
    Point to_local_point(CellCoordinate cc, CPosition pos) const {
        return {
            curve1.curve_length(pos[0]) - curve1.curve_length(cc[0]),
            curve2.curve_length(pos[1]) - curve2.curve_length(cc[1])
        };
    }

    /**
     * Convert a local cell point to a global CPosition.
     */
    CPosition to_cposition(CellCoordinate cc, Point p) const {
        return {{
            curve1.get_cpoint(cc[0], p.x),
            curve2.get_cpoint(cc[1], p.y)
        }};
    }

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

    template<ImageMetric imageMetric, ParamMetric paramMetric>
    Points compute_path(const Cell& cell, const Point& s, const Point& t) const;

    template<ImageMetric imageMetric, ParamMetric paramMetric>
    void steepest_descent(const Cell& cell, Point s, const Point& t, Points& path) const;

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

    Points compute_matching();

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(const Node& node, std::vector<Node>& neighbors) const;
    cost_t cost(const Node& s, const Node& t) const;
    cost_t heuristic_cost(const Node& s, const Node& goal) const;
};