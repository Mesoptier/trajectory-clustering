#pragma once

#include <Curve.h>
#include <geom.h>
#include <iostream>
#include <utility>
#include "metrics.h"

using CellCoordinate = std::array<PointID, 2>;

struct Cell {
    // Source: bottom-left corner
    const Point& s1;
    const Point& s2;

    // Target: top-right corner
    const Point& t1;
    const Point& t2;

    // Edge lengths
    distance_t len1;
    distance_t len2;

    // NOTE: These values are undefined for degenerate cells
    Point center;
    Line ell_m; // Monotone axis
//    Line ellH; // Line through points where ellipse tangent is horizontal
//    Line ellV; // Line through points where ellipse tangent is vertical

    Cell(const Point& s1, const Point& s2, const Point& t1, const Point& t2) :
        s1(s1), s2(s2), t1(t1), t2(t2),
        len1(s1.dist(t1)),
        len2(s2.dist(t2))
    {
        if (!is_degenerate()) {
            const Edge e1(s1, t1);
            const Edge e2(s2, t2);

            if (isParallel(e1.line, e2.line)) {
                const auto image_point = e1.line.closest(e2.first);
                center = {e1.param(image_point), 0};
            } else {
                const auto image_point = intersect(e1.line, e2.line);
                center = {e1.param(image_point), e2.param(image_point)};
            }

            ell_m = {center, {1, 1}};
        }
    }

    Point s() const {
        return {0, 0};
    }

    Point t() const {
        return {len1, len2};
    }

    bool is_degenerate() const {
        return is_e1_degenerate() || is_e2_degenerate();
    }
    bool is_e1_degenerate() const {
        return len1 == 0;
    }
    bool is_e2_degenerate() const {
        return len2 == 0;
    }

    /**
     * Get points in image space corresponding to the given point in
     * parameter space.
     *
     * @param p
     * @return
     */
    std::pair<Point, Point> interpolate_at(const Point& p) const {
        const Point p1 = is_e1_degenerate() ? s1 : Edge(s1, t1).interpolate_at(p.x);
        const Point p2 = is_e2_degenerate() ? s2 : Edge(s2, t2).interpolate_at(p.y);
        return {p1, p2};
    }
};

class IntegralFrechet
{
private:
    const Curve curve1;
    const Curve curve2;

    /**
     * Compute the cost of the optimal matching from the bottom-left corner to
     * the top-right corner of the cell.
     *
     * @tparam imageMetric
     * @tparam paramMetric
     * @param cell
     * @return
     */
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t cost(const Cell& cell) const;

    /**
     * Compute the optimal matching from the bottom-left corner to the top-right
     * corner of the cell.
     *
     * @tparam imageMetric
     * @tparam paramMetric
     * @param cell
     * @return
     */
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    Points compute_matching(const Cell& cell) const;
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    void steepest_descent(const Cell& cell, Point s, const Point& t, Points& path) const;

    // Compute cost over path by integration
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t integrate(const Cell& cell, const Points& cell_matching) const;
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t integrate(const Cell& cell, const Point& s, const Point& t) const;

    // TODO: Remove old integration methods
    template<ImageMetric imageMetric, ParamMetric paramMetric>
    distance_t integrate(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const;

    // Cost from distance between curves
    template<ImageMetric imageMetric>
    distance_t integrate_cost(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const;

    // Cost from distance travelled over curve
    template<ParamMetric paramMetric>
    distance_t integrate_dist(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const;

public:

    IntegralFrechet(const Curve& curve1, const Curve& curve2);

    std::pair<distance_t, Points> compute_matching();

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(const Node& node, std::vector<Node>& neighbors) const;
    cost_t cost(const Node& s, const Node& t) const;
    cost_t heuristic_cost(const Node& s, const Node& goal) const;

    Point node_as_point(const Node& node) const {
        return {
            curve1.curve_length(node[0]),
            curve2.curve_length(node[1])
        };
    }
};