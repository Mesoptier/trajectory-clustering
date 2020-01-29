#pragma once

#include <iostream>
#include <utility>
#include "metrics.h"
#include "../geom.h"
#include "../Curve.h"
#include "Cell.h"
#include "../a_star.h"

class IntegralFrechet
{
private:
    const Curve& curve1;
    const Curve& curve2;

    const ParamMetric param_metric;
    const distance_t  resolution;

    Cell get_cell(const CPosition& s, const CPosition& t) const;

    Points compute_cell_matching(const Cell& cell) const;
    distance_t compute_cell_cost(const Cell& cell, const Points& matching) const;

    /**
     * Compute the cost of the optimal matching from the bottom-left corner to
     * the top-right corner of the cell.
     *
     * @param cell
     * @return
     */
    distance_t cost(const Cell& cell) const;

public:

    IntegralFrechet(
        const Curve& curve1,
        const Curve& curve2,
        ParamMetric param_metric,
        distance_t resolution
    );

    struct MatchingResult {
        distance_t cost;
        Points matching;
        shortest_path_algs::SearchStat search_stat;
    };

    MatchingResult compute_matching();

    //
    // Requirements for A* algorithm:
    //

    using cost_t = distance_t;
    using Node = CPosition;
    void get_neighbors(const Node& node, std::vector<Node>& neighbors, BFDirection dir = BFDirection::Forward) const;
    cost_t cost(const Node& s, const Node& t) const;
    cost_t heuristic_cost(const Node& s, const Node& goal) const;

    Point node_as_point(const Node& node) const {
        return {
            curve1.curve_length(node[0]),
            curve2.curve_length(node[1])
        };
    }
};