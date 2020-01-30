#pragma once

#include "../Cell.h"

/**
 * Compute the optimal matching from the bottom-left corner to the top-right
 * corner of the cell.
 *
 * @tparam param_metric
 * @param cell
 * @return
 */
template<ParamMetric param_metric>
Points compute_matching(const Cell& cell, const Point& s, const Point& t);

/**
 * Compute the part of the cost from the distance between the curves.
 * Note: does not depend on the ParamMetric used.
 *
 * @param cell
 * @return
 */
distance_t integrate_linear_cost(const Cell& cell, const Point& s, const Point& t);

/**
 * Compute the part of the cost from distance travelled over the curve.
 *
 * @tparam param_metric
 * @param cell
 * @return
 */
template<ParamMetric param_metric>
distance_t integrate_linear_dist(const Cell& cell, const Point& s, const Point& t);

/**
 * Compute the cost of a single edge in parameter space.
 *
 * @tparam param_metric
 * @param cell
 * @return
 */
template<ParamMetric param_metric>
distance_t integrate_linear(const Cell& cell, const Point& s, const Point& t) {
    return integrate_linear_cost(cell, s, t) * integrate_linear_dist<param_metric>(cell, s, t);
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
        cost += integrate_linear<param_metric>(cell, matching[i - 1], matching[i]);
    }
    return cost;
}