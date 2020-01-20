#pragma once

#include "../geom.h"
#include "metrics.h"

namespace NewCell {
    struct Cell
    {
        // Bottom-left corner of the cell
        Point s1; // Image space: start of edge 1
        Point s2; // Image space: start of edge 2

        // Top-right corner of the cell
        Point t1; // Image space: end of edge 1
        Point t2; // Image space: end of edge 2

        distance_t len1;
        distance_t len2;

        Point s;
        Point t;

        // Parameters to the height function
        // Note: d(x,y) = (x-a)^2 + (y-b)^2 + l*(x-a)*(y-b) + c
        Point mid; // = (a, b)
        distance_t c;
        distance_t l;

        // Ellipse axes
        Line ell_m;
        Line ell_h;
        Line ell_v;

        Cell(const Point& s1, const Point& s2, const Point& t1, const Point& t2) :
            s1(s1), s2(s2), t1(t1), t2(t2),
            len1(s1.dist(t1)),
            len2(s2.dist(t2)),
            s({0, 0}),
            t({len1, len2})
        {
            if (!is_degenerate()) {
                const auto l1 = Line::fromTwoPoints(s1, t1);
                const auto l2 = Line::fromTwoPoints(s2, t2);

                // Find ellipse midpoint
                if (!isParallel(l1, l2)) {
                    const auto p = intersect(l1, l2);
                    mid = {l1(p), l2(p)};
                    c = 0;
                } else {
                    const auto p = l2.closest(s1);
                    mid = {0, l2(p)};
                    c = p.dist(s1);
                }

                // Compute ellipse axes
                ell_m = Line(mid, {1, 1});
                const auto v = dot(l2.direction, l1.direction);
                ell_h = Line(mid, {v, 1});
                ell_v = Line(mid, {1, v});
            }
        }

        bool is_degenerate() {
            return approx_equal(s1, t1) || approx_equal(s2, t2);
        }

        /**
         * Get points in image space corresponding to the given point in
         * parameter space.
         *
         * @param p
         * @return
         */
        [[nodiscard]]
        std::pair<Point, Point> interpolate_at(const Point& p) const;

        /**
         * Get a subcell of this cell.
         *
         * @param s Bottom-left corner (in coordinates relative to this cell) of the subcell.
         * @param t Top-right corner (in coordinates relative to this cell) of the subcell.
         * @return
         */
        [[nodiscard]]
        Cell subcell(const Point& s, const Point& t) const;
    };

    // TODO: Move methods out of Cell namespace

    /**
     * Compute the optimal matching from the bottom-left corner to the top-right
     * corner of the cell.
     *
     * @tparam param_metric
     * @param cell
     * @return
     */
    template<ParamMetric param_metric>
    Points compute_matching(const Cell& cell);

    template<ParamMetric param_metric>
    Points steepest_descent(const Cell& cell, const Point& s, const Point& t, BFDirection dir);

    /**
     * Compute the part of the cost from the distance between the curves.
     * Note: does not depend on the ParamMetric used.
     *
     * @param cell
     * @return
     */
    distance_t integrate_linear_cost(const Cell& cell);

    /**
     * Compute the part of the cost from distance travelled over the curve.
     *
     * @tparam param_metric
     * @param cell
     * @return
     */
    template<ParamMetric param_metric>
    distance_t integrate_linear_dist(const Cell& cell);

    /**
     * Compute the cost of a single edge in parameter space.
     *
     * @tparam param_metric
     * @param cell
     * @return
     */
    template<ParamMetric param_metric>
    distance_t integrate_linear(const Cell& cell) {
        return integrate_linear_cost(cell) * integrate_linear_dist<param_metric>(cell);
    }

    /**
     * Compute the cost of the matching in the cell.
     *
     * @tparam param_metric
     * @param cell
     * @param matching
     * @return
     */
    template<ParamMetric param_metric>
    distance_t compute_cost(const Cell& cell, const Points& matching) {
        distance_t cost = 0;
        for (size_t i = 1; i < matching.size(); ++i) {
            cost += NewCell::integrate_linear<param_metric>(
                cell.subcell(matching[i - 1], matching[i])
            );
        }
        return cost;
    }
}
