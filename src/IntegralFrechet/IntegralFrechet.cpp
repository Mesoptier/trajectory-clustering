#include <utility>
#include "IntegralFrechet.h"
#include "a_star.h"

IntegralFrechet::IntegralFrechet(Curve curve1, Curve curve2) : curve1(std::move(curve1)), curve2(std::move(curve2)) {}

void IntegralFrechet::findPath() {
    Node start{CPoint(0, 0), CPoint(0, 0)};
    Node goal{CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)};

    a_star_search(*this, start, goal);
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
    return 0;
}

IntegralFrechet::cost_t
IntegralFrechet::heuristic_cost(const IntegralFrechet::Node& s, const IntegralFrechet::Node& goal) const {
    return 0;
}
