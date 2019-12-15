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

Points IntegralFrechet::compute_matching() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    auto cmatching = a_star_search(*this, start, goal);

    Points matching; // Arc-length parametrized
    for (auto cpos : cmatching) {
        matching.emplace_back(curve1.curve_length(cpos[0]), curve2.curve_length(cpos[1]));
    }
    return matching;
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
Points IntegralFrechet::compute_path<ImageMetric::L2_Squared, ParamMetric::L1>(const Cell& cell, const Point& s, const Point& t) const {


    // TODO: Compute actual optimal path
    return {s, t};
}

template<>
CPositions IntegralFrechet::steepest_descent<ImageMetric::L2_Squared, ParamMetric::L1>(const CPosition& s, const CPosition& t) const {
    return CPositions();
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