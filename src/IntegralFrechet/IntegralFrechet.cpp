#include "IntegralFrechet/IntegralFrechet.h"

#include "utils/a_star.h"
#include "IntegralFrechet/Cell.h"
#include "IntegralFrechet/metrics/include.h"

IntegralFrechet::IntegralFrechet(
    Curve const& c1,
    Curve const& c2,
    ParamMetric metric,
    distance_t res,
    MatchingBand const* const b
) :
    curve1(c1),
    curve2(c2),
    resolution(res),
    band(b),
    param_metric(metric) {}

Cell IntegralFrechet::get_cell(CPosition const& s, CPosition const& t) const {
    auto const s1 = curve1.interpolate_at(s[0]);
    auto const s2 = curve2.interpolate_at(s[1]);
    auto const t1 = curve1.interpolate_at(t[0]);
    auto const t2 = curve2.interpolate_at(t[1]);
    return {s1, s2, t1, t2};
}

IntegralFrechet::MatchingResult IntegralFrechet::compute_matching() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    // Matching from nodes to nodes (on cell boundaries)
    auto const search_result = bidirectional_dijkstra_search(*this, start, goal);
    auto const& node_matching = search_result.path;

    // Actual polygonal matching with arc-length coordinates
    Points matching{{0, 0}};

    for (std::size_t i = 1; i < node_matching.size(); ++i) {
        auto const s = node_matching[i - 1];
        auto const t = node_matching[i];

        auto const cell = get_cell(s, t);
        auto cell_matching = compute_cell_matching(cell, cell.s, cell.t);

        Point const offset(
            curve1.curve_length(s[0]),
            curve2.curve_length(s[1])
        );

        for (std::size_t j = 1; j < cell_matching.size(); ++j)
            matching.push_back(offset + cell_matching[j]);
    }
    return {
        search_result.cost,
        matching,
        search_result.stat,
    };
}

distance_t IntegralFrechet::cost(Cell const& cell,
        Point const& s, Point const& t) const {
    Points const cell_matching = compute_cell_matching(cell, s, t);
    return compute_cell_cost(cell, cell_matching);
}

extern template
Points compute_matching<ParamMetric::L1>(Cell const&, Point const&, Point const&);
extern template
Points compute_matching<ParamMetric::LInfinity_NoShortcuts>(Cell const&,
    Point const&, Point const&);

Points IntegralFrechet::compute_cell_matching(Cell const& cell,
        Point const& s, Point const& t) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_matching<ParamMetric::L1>(cell, s, t);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_matching<ParamMetric::LInfinity_NoShortcuts>(
                cell, s, t);
        default:
            throw std::invalid_argument("Unsupported norm");
    }
}

extern template
distance_t integrate_linear_dist<ParamMetric::L1>(Cell const&, Point const& s,
    Point const& t);
extern template
distance_t integrate_linear_dist<ParamMetric::LInfinity_NoShortcuts>(
    Cell const&, Point const& s, Point const& t);

distance_t IntegralFrechet::compute_cell_cost(Cell const& cell,
        Points const& matching) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_cost<ParamMetric::L1>(cell, matching);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_cost<ParamMetric::LInfinity_NoShortcuts>(
                cell, matching);
        default:
            throw std::invalid_argument("Unsupported norm");
    }
}

//
// Requirements for A* algorithm:
//
void IntegralFrechet::get_neighbors(Node const& node,
        std::vector<Node>& neighbors, std::vector<cost_t>& costs,
        BFDirection dir) const {
    // p: point in the cell
    CPoint const p1 = node[0];
    CPoint const p2 = node[1];

    // s: corner of the cell smaller than p w.r.t. `dir` (i.e. s <=_dir p)
    CPoint const s1 = dir == BFDirection::Forward ? p1.floor() : p1.ceil();
    CPoint const s2 = dir == BFDirection::Forward ? p2.floor() : p2.ceil();

    // t: corner of the cell greater than p w.r.t. `dir` (i.e. p <=_dir t)
    CPoint const t1 = dir == BFDirection::Forward
        ? (p1.getPoint() == curve1.size() - 1 ? p1 : p1.floor() + 1)
        : (p1.ceil().getPoint() == 0 ? p1 : p1.ceil() - 1);
    CPoint const t2 = dir == BFDirection::Forward
        ? (p2.getPoint() == curve2.size() - 1 ? p2 : p2.floor() + 1)
        : (p2.ceil().getPoint() == 0 ? p2 : p2.ceil() - 1);

    //
    // COMPUTE NEIGHBORS
    //

    if (p1 == t1 && p2 == t2) {
        // Subcell has no width and no height (Degenerate case #1)
        // -> We are at the end of the parameter region, nowhere to go from here.
        return;
    }

    // To prevent floating point errors we need to have a full cell in the forwards direction
    // and compute sampling points along its boundary.
    Cell const full_cell = dir == BFDirection::Forward
        ? get_cell({s1, s2}, {t1, t2})
        : get_cell({t1, t2}, {s1, s2});

    // CPosition of bottom-left corner of full_cell
    CPoint const fs1 = dir == BFDirection::Forward ? s1 : t1;
    CPoint const fs2 = dir == BFDirection::Forward ? s2 : t2;

    if (p1 == t1 || p2 == t2) {
        // Subcell has no width OR no height (Degenerate case #2)
        // -> We are along the parameter region boundary, only one neighbor makes sense.
        if (band == nullptr || band->contains_point(t1, t2))
            neighbors.push_back({t1, t2});
    }
    else {
        // At this point we have s <=_{dir} p <_{dir} t

        // Number of nodes along horizontal/vertical boundaries
        auto const n_h = static_cast<std::size_t>(
            std::ceil(full_cell.len1 / resolution)) + 1;
        auto const n_v = static_cast<std::size_t>(
            std::ceil(full_cell.len2 / resolution)) + 1;

        // TODO: Add intersections of ell_h/ell_v/ell_m with cell boundaries
        // (and of adjacent cells)
        // Cell boundary lines
        // Line const bound_h = Line::horizontal(dir == BFDirection::Forward
        //     ? full_cell.t : full_cell.s);
        // Line const bound_v = Line::vertical(dir == BFDirection::Forward
        //     ? full_cell.t : full_cell.s);

        // If there another cell above/below this cell...
        if (dir == BFDirection::Forward ? (t2.getPoint() + 1 < curve2.size())
                : (t2.getPoint() >= 1)) {
            // ...sample points along the top/bottom boundary of the cell:
            for (std::size_t i = 0; i < n_h; ++i) {
                CPoint const o1 = fs1 + (static_cast<distance_t>(i) / (n_h - 1.0));
                if (dir == BFDirection::Forward ? o1 < p1 : o1 > p1)
                    continue;
                if (band != nullptr && !band->contains_point(o1, t2))
                    continue;

                neighbors.push_back({o1, t2});
            }
        }

        // If there another cell left/right of this cell...
        if (dir == BFDirection::Forward ? (t1.getPoint() + 1 < curve1.size())
                : (t1.getPoint() >= 1)) {
            // ...sample points along the left/right boundary of the cell:
            for (std::size_t i = 0; i < n_v; ++i) {
                CPoint const o2 = fs2 + (static_cast<distance_t>(i) / (n_v - 1.0));
                if (dir == BFDirection::Forward ? o2 < p2 : o2 > p2)
                    continue;
                if (band != nullptr && !band->contains_point(t1, o2))
                    continue;

                neighbors.push_back({t1, o2});
            }
        }

        // Always try to add top-right (or bottom-left) corner
        if (band == nullptr || band->contains_point(t1, t2))
            neighbors.push_back({t1, t2});
    }

    //
    // COMPUTE COSTS
    //
    Point const s = {
        curve1.curve_length(fs1, node[0]),
        curve2.curve_length(fs2, node[1]),
    };

    for (auto const& neighbor : neighbors) {
        Point const t = {
            curve1.curve_length(fs1, neighbor[0]),
            curve2.curve_length(fs2, neighbor[1]),
        };

        if (dir == BFDirection::Forward)
            costs.push_back(cost(full_cell, s, t));
        else
            costs.push_back(cost(full_cell, t, s));
    }
}
