#include "IntegralFrechet.h"

#include "a_star.h"
#include "Cell.h"
#include "metrics/include.h"

IntegralFrechet::IntegralFrechet(
    const Curve& curve1,
    const Curve& curve2,
    ParamMetric param_metric,
    distance_t resolution
) : curve1(curve1), curve2(curve2), param_metric(param_metric), resolution(resolution) {}

Cell IntegralFrechet::get_cell(const CPosition& s, const CPosition& t) const {
    const auto s1 = curve1.interpolate_at(s[0]);
    const auto s2 = curve2.interpolate_at(s[1]);
    const auto t1 = curve1.interpolate_at(t[0]);
    const auto t2 = curve2.interpolate_at(t[1]);
    return {s1, s2, t1, t2};
}

std::pair<distance_t, Points> IntegralFrechet::compute_matching() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    // Matching from nodes to nodes (on cell boundaries)
    std::vector<Node> node_matching;
    cost_t cost;
    std::tie(cost, node_matching) = a_star_search(*this, start, goal);

    // Actual polygonal matching with arc-length coordinates
    Points matching{{0, 0}};

    for (int i = 1; i < node_matching.size(); ++i) {
        const auto s = node_matching[i - 1];
        const auto t = node_matching[i];

        const auto cell = get_cell(s, t);
        auto cell_matching = compute_cell_matching(cell);

        const Point offset(
            curve1.curve_length(s[0]),
            curve2.curve_length(s[1])
        );

        for (size_t j = 1; j < cell_matching.size(); ++j) {
            matching.push_back(offset + cell_matching[j]);
        }
    }
    return {cost, matching};
}

distance_t IntegralFrechet::cost(const Cell& cell) const {
    const Points cell_matching = compute_cell_matching(cell);
    return compute_cell_cost(cell, cell_matching);
}

Points IntegralFrechet::compute_cell_matching(const Cell& cell) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_matching<ParamMetric::L1>(cell);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell);
    }
}

distance_t IntegralFrechet::compute_cell_cost(const Cell& cell, const Points& matching) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_cost<ParamMetric::L1>(cell, matching);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_cost<ParamMetric::LInfinity_NoShortcuts>(cell, matching);
    }
}

//
// Requirements for A* algorithm:
//

void IntegralFrechet::get_neighbors(const IntegralFrechet::Node& node, std::vector<Node>& neighbors) const {
    if (node[0].getPoint() + 1 == curve1.size()) { // Right boundary
        if (node[1].getPoint() + 1 == curve2.size()) { // Top-right corner
            return; // No neighbors
        }
        // Move up one cell
        neighbors.push_back({node[0], node[1].floor() + 1});
        return;
    }
    if (node[1].getPoint() + 1 == curve2.size()) { // Top boundary
        // Move right one cell
        neighbors.push_back({node[0].floor() + 1, node[1]});
        return;
    }

    // We are not along the outer top/right boundary, which means the node is part of a cell

    // Compute the number of sampling points
    const auto len1 = curve1.curve_length(node[0], node[0].floor() + 1);
    const auto len2 = curve2.curve_length(node[1], node[1].floor() + 1);
    const size_t n1 = ceil(len1 / resolution) + 1;
    const size_t n2 = ceil(len2 / resolution) + 1;

    // Sample along top edge if there are more cells above this one
    if (node[1].getPoint() + 2 < curve2.size()) {
        for (int i = 0; i < n1 - 1; ++i) {
            distance_t fraction = i / (n1 - 1.);
            if (fraction >= node[0].getFraction()) {
                neighbors.push_back({node[0].floor() + fraction, node[1].floor() + 1});
            }
        }
    }

    // Sample along right edge if there are more cells right of this one
    if (node[0].getPoint() + 2 < curve1.size()) {
        for (int i = 0; i < n2 - 1; ++i) {
            distance_t fraction = i / (n2 - 1.);
            if (fraction >= node[1].getFraction()) {
                neighbors.push_back({node[0].floor() + 1, node[1].floor() + fraction});
            }
        }
    }

    // Always add top-right corner
    neighbors.push_back({node[0].floor() + 1, node[1].floor() + 1});
}

IntegralFrechet::cost_t IntegralFrechet::cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& t) const {
    return cost(get_cell(s, t));
}

IntegralFrechet::cost_t
IntegralFrechet::heuristic_cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& goal) const {
    return 0;
}