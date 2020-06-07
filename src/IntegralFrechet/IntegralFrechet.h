#pragma once

#include <iostream>
#include <utility>
#include "IntegralFrechet/metrics.h"
#include "geom.h"
#include "Curve.h"
#include "IntegralFrechet/Cell.h"
#include "utils/a_star.h"
#include "IntegralFrechet/MatchingBand.h"

class IntegralFrechet
{
private:
    const Curve& curve1;
    const Curve& curve2;
    const distance_t  resolution;
    const MatchingBand* const band;
    const ParamMetric param_metric;

    Cell get_cell(const CPosition& s, const CPosition& t) const;

    Points compute_cell_matching(const Cell& cell, const Point& s, const Point& t) const;
    distance_t compute_cell_cost(const Cell& cell, const Points& matching) const;

    /**
     * Compute the cost of the optimal matching from the bottom-left corner to
     * the top-right corner of the cell.
     */
    distance_t cost(const Cell& cell, const Point& s, const Point& t) const;

public:

    IntegralFrechet(const Curve& c1, const Curve& c2, ParamMetric metric,
        distance_t res, const MatchingBand* b = nullptr);

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
    void get_neighbors(const Node& node, std::vector<Node>& neighbors, std::vector<cost_t>& costs, BFDirection dir) const;

    Point node_as_point(const Node& node) const {
        return {
            curve1.curve_length(node[0]),
            curve2.curve_length(node[1])
        };
    }
};
