#include "IntegralFrechet.h"
#include "a_star.h"

IntegralFrechet::IntegralFrechet(const Curve& curve1, const Curve& curve2) : curve1(curve1), curve2(curve2) {
    cells.reserve(curve1.size() * curve2.size());

    for (PointID p1 = 0; p1 + 1 < curve1.size(); ++p1) {
        for (PointID p2 = 0; p2 + 1 < curve2.size(); ++p2) {
            // TODO: Proper initialization of cells
            cells.emplace_back();
        }
    }
}

CPositions IntegralFrechet::compute_matching() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    return a_star_search(*this, start, goal);
}

const Cell& IntegralFrechet::get_cell(PointID p1, PointID p2) const {
    return cells[p1 * curve2.size() + p2];
}

template<ImageMetric imageMetric, ParamMetric paramMetric>
distance_t IntegralFrechet::cost(const CPosition& s, const CPosition& t) const {
    auto path = compute_path<imageMetric, paramMetric>(s, t);
    return integrate<imageMetric, paramMetric>(path);
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
CPositions IntegralFrechet::compute_path<ImageMetric::L2_Squared, ParamMetric::L1>(const CPosition& s, const CPosition& t) const {
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
    if (node[0] < curve1.size() - 1) {
        neighbors.push_back({node[0] + 1, node[1]});

        if (node[1] < curve2.size() - 1) {
            neighbors.push_back({node[0] + 1, node[1] + 1});
        }
    }
    if (node[1] < curve2.size() - 1) {
        neighbors.push_back({node[0], node[1] + 1});
    }
}

IntegralFrechet::cost_t IntegralFrechet::cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& t) const {
    return cost<ImageMetric::L2_Squared, ParamMetric::L1>(s, t);
}

IntegralFrechet::cost_t
IntegralFrechet::heuristic_cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& goal) const {
    return 0;
}