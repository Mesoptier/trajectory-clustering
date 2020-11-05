#ifndef INTEGRAL_FRECHET_H
#define INTEGRAL_FRECHET_H

#include <iostream>
#include <utility>
#include "IntegralFrechet/metrics.h"
#include "geom.h"
#include "Curve.h"
#include "IntegralFrechet/Cell.h"
#include "utils/a_star.h"
#include "IntegralFrechet/MatchingBand.h"

class IntegralFrechet {
private:
    Curve const& curve1;
    Curve const& curve2;
    distance_t const resolution;
    MatchingBand const* const band;
    ParamMetric const param_metric;

    Cell get_cell(CPosition const& s, CPosition const& t) const;

    Points compute_cell_matching(Cell const& cell,
        Point const& s, Point const& t) const;
    distance_t compute_cell_cost(Cell const& cell,
        Points const& matching) const;

    /**
     * Compute the cost of the optimal matching from the bottom-left corner to
     * the top-right corner of the cell.
     */
    distance_t cost(Cell const& cell, Point const& s, Point const& t) const;

public:

    IntegralFrechet(Curve const& c1, Curve const& c2, ParamMetric metric,
        distance_t res, MatchingBand const* b = nullptr);

    struct MatchingResult {
        distance_t cost;
        Points matching;
        SearchStat search_stat;
    };

    MatchingResult compute_matching();

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(Node const& node, std::vector<Node>& neighbors,
        std::vector<cost_t>& costs, BFDirection dir) const;

    Point node_as_point(Node const& node) const {
        return {
            curve1.curve_length(node[0]),
            curve2.curve_length(node[1])
        };
    }
};
#endif
