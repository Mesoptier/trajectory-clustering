#include "IntegralFrechet.h"

#include "a_star.h"
#include "Cell.h"
#include "metrics/include.h"

IntegralFrechet::IntegralFrechet(
    const Curve& c1,
    const Curve& c2,
    ParamMetric metric,
    distance_t res,
    const MatchingBand* const b
) :
    curve1(c1),
    curve2(c2),
    resolution(res),
    band(b),
    param_metric(metric) {}

Cell IntegralFrechet::get_cell(const CPosition& s, const CPosition& t) const {
    const auto s1 = curve1.interpolate_at(s[0]);
    const auto s2 = curve2.interpolate_at(s[1]);
    const auto t1 = curve1.interpolate_at(t[0]);
    const auto t2 = curve2.interpolate_at(t[1]);
    // std::cout << s1 << "\n" << s2 << "\n" << t1 << "\n" << t2 << "\n";
    // std::cout << s[0] << "\n" << t[0] << "\n";
    // std::cout << s1 << "\n" << t1 << "\n"; 
    // std::cout << s1.dist_sqr(t1) << "\n";
    return {s1, s2, t1, t2};
}

IntegralFrechet::MatchingResult IntegralFrechet::compute_matching() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    // Matching from nodes to nodes (on cell boundaries)
    const auto search_result = bidirectional_dijkstra_search(*this, start, goal);
    const auto& node_matching = search_result.path;

    // Actual polygonal matching with arc-length coordinates
    Points matching{{0, 0}};

    for (size_t i = 1; i < node_matching.size(); ++i) {
        const auto s = node_matching[i - 1];
        const auto t = node_matching[i];

        const auto cell = get_cell(s, t);
        auto cell_matching = compute_cell_matching(cell, cell.s, cell.t);

        const Point offset(
            curve1.curve_length(s[0]),
            curve2.curve_length(s[1])
        );

        for (size_t j = 1; j < cell_matching.size(); ++j) {
            matching.push_back(offset + cell_matching[j]);
        }
    }
    return {
        search_result.cost,
        matching,
        search_result.stat,
    };
}

distance_t IntegralFrechet::cost(const Cell& cell, const Point& s, const Point& t) const {
    const Points cell_matching = compute_cell_matching(cell, s, t);
    return compute_cell_cost(cell, cell_matching);
}

extern template
Points compute_matching<ParamMetric::L1>(const Cell&, const Point&, const Point&);
extern template
Points compute_matching<ParamMetric::LInfinity_NoShortcuts>(const Cell&,
    const Point&, const Point&);

Points IntegralFrechet::compute_cell_matching(const Cell& cell, const Point& s, const Point& t) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_matching<ParamMetric::L1>(cell, s, t);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell, s, t);
        default:
            throw std::invalid_argument("Unsupported norm");
    }
}

extern template
distance_t integrate_linear_dist<ParamMetric::L1>(const Cell&, const Point& s,
    const Point& t);
extern template
distance_t integrate_linear_dist<ParamMetric::LInfinity_NoShortcuts>(
    const Cell&, const Point& s, const Point& t);

distance_t IntegralFrechet::compute_cell_cost(const Cell& cell, const Points& matching) const {
    switch (param_metric) {
        case ParamMetric::L1:
            return ::compute_cost<ParamMetric::L1>(cell, matching);
        case ParamMetric::LInfinity_NoShortcuts:
            return ::compute_cost<ParamMetric::LInfinity_NoShortcuts>(cell, matching);
        default:
            throw std::invalid_argument("Unsupported norm");
    }
}

//
// Requirements for A* algorithm:
//

void IntegralFrechet::get_neighbors(const Node& node, std::vector<Node>& neighbors, std::vector<cost_t>& costs, BFDirection dir) const {
    // p: point in the cell
    const CPoint p1 = node[0];
    const CPoint p2 = node[1];

    // s: corner of the cell smaller than p w.r.t. `dir` (i.e. s <=_dir p)
    const CPoint s1 = dir == BFDirection::Forward ? p1.floor() : p1.ceil();
    const CPoint s2 = dir == BFDirection::Forward ? p2.floor() : p2.ceil();

    // t: corner of the cell greater than p w.r.t. `dir` (i.e. p <=_dir t)
    const CPoint t1 = dir == BFDirection::Forward
        ? (p1.getPoint() == curve1.size() - 1 ? p1 : p1.floor() + 1)
        : (p1.ceil().getPoint() == 0 ? p1 : p1.ceil() - 1);
    const CPoint t2 = dir == BFDirection::Forward
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
    const Cell full_cell = dir == BFDirection::Forward
        ? get_cell({s1, s2}, {t1, t2})
        : get_cell({t1, t2}, {s1, s2});

    // CPosition of bottom-left corner of full_cell
    const CPoint fs1 = dir == BFDirection::Forward ? s1 : t1;
    const CPoint fs2 = dir == BFDirection::Forward ? s2 : t2;

    if (p1 == t1 || p2 == t2) {
        // Subcell has no width OR no height (Degenerate case #2)
        // -> We are along the parameter region boundary, only one neighbor makes sense.
        if (band == nullptr || band->contains_point(t1, t2)) {
            neighbors.push_back({t1, t2});
        }
    } else {
        // At this point we have s <=_{dir} p <_{dir} t

        // Number of nodes along horizontal/vertical boundaries
        const size_t n_h = static_cast<size_t>(std::ceil(full_cell.len1 / resolution)) + 1;
        const size_t n_v = static_cast<size_t>(std::ceil(full_cell.len2 / resolution)) + 1;

        if (n_h > 10 || n_v  > 10) {
            std::cout << "this is weird...\n";
        }

        // TODO: Add intersections of ell_h/ell_v/ell_m with cell boundaries (and of adjacent cells)
        // Cell boundary lines
        // const Line bound_h = Line::horizontal(dir == BFDirection::Forward ? full_cell.t : full_cell.s);
        // const Line bound_v = Line::vertical(dir == BFDirection::Forward ? full_cell.t : full_cell.s);

        // If there another cell above/below this cell...
        if (dir == BFDirection::Forward ? (t2.getPoint() + 1 < curve2.size()) : (t2.getPoint() >= 1)) {
            // ...sample points along the top/bottom boundary of the cell:
            for (size_t i = 0; i < n_h; ++i) {
                const CPoint o1 = fs1 + (static_cast<distance_t>(i) / (n_h - 1.0));
                if (dir == BFDirection::Forward ? o1 < p1 : o1 > p1) {
                    continue;
                }
                if (band != nullptr && !band->contains_point(o1, t2)) {
                    continue;
                }

                neighbors.push_back({o1, t2});
            }
        }

        // If there another cell left/right of this cell...
        if (dir == BFDirection::Forward ? (t1.getPoint() + 1 < curve1.size()) : (t1.getPoint() >= 1)) {
            // ...sample points along the left/right boundary of the cell:
            for (size_t i = 0; i < n_v; ++i) {
                const CPoint o2 = fs2 + (static_cast<distance_t>(i) / (n_v - 1.0));
                if (dir == BFDirection::Forward ? o2 < p2 : o2 > p2) {
                    continue;
                }
                if (band != nullptr && !band->contains_point(t1, o2)) {
                    continue;
                }

                neighbors.push_back({t1, o2});
            }
        }

        // Always try to add top-right (or bottom-left) corner
        if (band == nullptr || band->contains_point(t1, t2)) {
            neighbors.push_back({t1, t2});
        }
    }

    //
    // COMPUTE COSTS
    //

    const Point s = {
        curve1.curve_length(fs1, node[0]),
        curve2.curve_length(fs2, node[1]),
    };

    for (const auto& neighbor : neighbors) {
        const Point t = {
            curve1.curve_length(fs1, neighbor[0]),
            curve2.curve_length(fs2, neighbor[1]),
        };

        if (dir == BFDirection::Forward) {
            costs.push_back(cost(full_cell, s, t));
        } else {
            costs.push_back(cost(full_cell, t, s));
        }
    }
}
