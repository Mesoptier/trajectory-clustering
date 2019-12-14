#include <utility>
#include "IntegralFrechet.h"
#include "a_star.h"

IntegralFrechet::IntegralFrechet(Curve curve1, Curve curve2) : curve1(std::move(curve1)), curve2(std::move(curve2)) {}

void IntegralFrechet::findPath() {
    Graph graph(curve1, curve2);
    Graph::Node start{0, 0};
    Graph::Node goal{curve1.size() - 1, curve2.size() - 1};

    a_star_search<Graph>(graph, start, goal);
}