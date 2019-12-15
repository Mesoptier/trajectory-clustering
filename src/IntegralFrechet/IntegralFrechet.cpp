#include "IntegralFrechet.h"
#include "a_star.h"

IntegralFrechet::IntegralFrechet(const Curve& curve1, const Curve& curve2) : curve1(curve1), curve2(curve2) {
    cells.reserve(curve1.size() * curve2.size());

    for (PointID p2 = 0; p2 + 1 < curve2.size(); ++p2) {
        for (PointID p1 = 0; p1 + 1 < curve1.size(); ++p1) {
            cells.emplace_back(curve1[p1], curve2[p2], curve1[p1 + 1], curve2[p2 + 1]);
        }
    }
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
        auto cell_matching = compute_path<ImageMetric::L2_Squared, ParamMetric::L1>(node_matching[i - 1], node_matching[i]);

        for (int j = 1; j < cell_matching.size(); ++j) {
            matching.emplace_back(
                curve1.curve_length(cell_matching[j][0]),
                curve2.curve_length(cell_matching[j][1])
            );
        }
    }
    return {cost, matching};
}

const Cell& IntegralFrechet::get_cell(CellCoordinate cc) const {
    assert(cc[0] + 1 < curve1.size() && cc[1] + 1 < curve2.size());
    return cells[cc[1] * (curve1.size() - 1) + cc[0]];
}

template<ImageMetric imageMetric, ParamMetric paramMetric>
distance_t IntegralFrechet::cost(const CPosition& s, const CPosition& t) const {
    auto path = compute_path<imageMetric, paramMetric>(s, t);
    return integrate<imageMetric, paramMetric>(path);
}

template<ImageMetric imageMetric, ParamMetric paramMetric>
CPositions IntegralFrechet::compute_path(const CPosition& s, const CPosition& t) const {
    if (s[0] == t[0] || s[1] == t[1]) {
        // s and t lie on the same horizontal/vertical axis, so there is only 1 possible path
        return { s, t };
    }

    const CellCoordinate cc{s[0].getPoint(), s[1].getPoint()};
    const auto cell = get_cell(cc);

    // Convert to local cell coordinates
    const Point sl = to_local_point(cc, s);
    const Point tl = to_local_point(cc, t);

    // Compute optimal path through cell
    const Points local_path = compute_path<imageMetric, paramMetric>(cell, sl, tl);

    // Convert path back to global coordinates
    CPositions path;
    path.reserve(local_path.size());
    for (const auto& p : local_path) {
        path.push_back(to_cposition(cc, p));
    }

    return path;
}

template<ImageMetric imageMetric, ParamMetric paramMetric>
Points IntegralFrechet::compute_path(const Cell& cell, const Point& s, const Point& t) const {
    Points path1;
    Points path2;

    steepest_descent<imageMetric, paramMetric>(cell, s, t, path1);
    steepest_descent<imageMetric, paramMetric>(cell, t, s, path2);

    assert(approx_equal(path1.back(), path2.back()));

    // Combine the two steepest descent paths into one, skipping the duplicated minimum
    path1.insert(path1.end(), path2.rbegin() + 1, path2.rend());
    return path1;
}

template<ImageMetric imageMetric, ParamMetric paramMetric>
distance_t IntegralFrechet::integrate(const CPositions& path) const {
    distance_t cost = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        auto s1 = curve1.interpolate_at(path[i - 1][0]);
        auto s2 = curve2.interpolate_at(path[i - 1][1]);
        auto t1 = curve1.interpolate_at(path[i][0]);
        auto t2 = curve2.interpolate_at(path[i][1]);

        cost += integrate_cost<imageMetric>(s1, s2, t1, t2) * integrate_dist<paramMetric>(s1, s2, t1, t2);
    }
    return cost;
}

template<>
distance_t IntegralFrechet::integrate_cost<ImageMetric::L2_Squared>(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const {
    // Get difference vectors between start- and endpoints of the two sub-edges.
    const auto d1 = s1 - s2;
    const auto d2 = t1 - t2;

    const distance_t dx1 = d1.x;
    const distance_t dy1 = d1.y;
    const distance_t dx2 = d2.x;
    const distance_t dy2 = d2.y;

    const distance_t a = pow(dx1 - dx2, 2) + pow(dy1 - dy2, 2);
    const distance_t b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    const distance_t c = dx1 * dx1 + dy1 * dy1;

    return (a / 3 + b / 2 + c);
}

template<>
distance_t IntegralFrechet::integrate_dist<ParamMetric::L1>(const Point& s1, const Point& s2, const Point& t1, const Point& t2) const {
    return s1.dist(t1) + s2.dist(t2);
}

//
// L2_Squared + L1
//

template<>
void IntegralFrechet::steepest_descent<ImageMetric::L2_Squared, ParamMetric::L1>(const Cell& cell, Point s, const Point& t, Points& path) const {
    path.push_back(s);

    // If source is monotone greater than target, we need to do steepest
    // descent in reverse monotone direction.
    auto dir = getMonotoneDirection(s, t);
    MonotoneComparator compare(dir);

    // Which return value of Line::side indicates that a point is on the left
    int line_left = (dir == BFDirection::Forward ? -1 : 1);

    auto ell_m_side = cell.ell_m.side(s);

    if (ell_m_side != 0) { // -> s is not on ell_m
        if (ell_m_side == line_left) { // -> s is left of ell_m
            s = std::min({
                intersect(Line::horizontal(s), cell.ell_m),
                intersect(Line::horizontal(s), Line::vertical(t))
            }, compare);
        } else  { // -> s is right of ell_m
            s = std::min({
                intersect(Line::vertical(s), cell.ell_m),
                intersect(Line::vertical(s), Line::horizontal(t))
            }, compare);
        }

        if (approx_equal(path.back(), s))
            return;
        path.push_back(s);

        ell_m_side = cell.ell_m.side(s);
    }

    if (ell_m_side == 0) { // -> s is on ell_m
        if (compare(s, cell.center)) {
            // Move along ell_m until center or until the point where steepest descent from t hits ell_m
            s = std::min({
                cell.center,
                intersect(cell.ell_m, Line::vertical(t)),
                intersect(cell.ell_m, Line::horizontal(t)),
            }, compare);

            if (approx_equal(path.back(), s))
                return;
            path.push_back(s);
        }
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

    CellCoordinate cc{node[0].getPoint(), node[1].getPoint()};
    auto cell = get_cell(cc);

    // Sample along top edge if there are more cells above this one
    if (node[1].getPoint() + 2 < curve2.size()) {
        for (int i = 0; i < cell.n1 - 1; ++i) {
            distance_t fraction = i / (cell.n1 - 1.);
            if (fraction >= node[0].getFraction()) {
                neighbors.push_back({node[0].floor() + fraction, node[1].floor() + 1});
            }
        }
    }

    // Sample along right edge if there are more cells right of this one
    if (node[0].getPoint() + 2 < curve1.size()) {
        for (int i = 0; i < cell.n2 - 1; ++i) {
            distance_t fraction = i / (cell.n2 - 1.);
            if (fraction >= node[1].getFraction()) {
                neighbors.push_back({node[0].floor() + 1, node[1].floor() + fraction});
            }
        }
    }

    // Always add top-right corner
    neighbors.push_back({node[0].floor() + 1, node[1].floor() + 1});
}

IntegralFrechet::cost_t IntegralFrechet::cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& t) const {
    return cost<ImageMetric::L2_Squared, ParamMetric::L1>(s, t);
}

IntegralFrechet::cost_t
IntegralFrechet::heuristic_cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& goal) const {
    return 0;
}