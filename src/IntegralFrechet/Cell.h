#pragma once

#include "../geom.h"
#include "metrics.h"

namespace NewCell {
    struct Cell
    {
        // Bottom-left corner of the cell
        const Point s1;
        const Point s2;

        // Top-right corner of the cell
        const Point t1;
        const Point t2;

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
