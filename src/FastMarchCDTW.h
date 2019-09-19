#ifndef CODE_FASTMARCHCDTW_H
#define CODE_FASTMARCHCDTW_H

#include <utility>
#include <vector>
#include "Curve.h"

namespace FastMarchCDTW {

    typedef std::pair<unsigned int, unsigned int> Point;

    struct Node {
        Node(Point point, double cost) : point(std::move(point)), cost(cost) {}

        Point point;
        double cost;
    };

    struct CompareNode {
        bool operator()(const Node& a, const Node& b) const {
            return a.cost > b.cost;
        }
    };

    enum Tag {
        Accepted,
        Considered,
        Far,
    };

    void compute(const Curve<double>& curve1, const Curve<double>& curve2);

}


#endif //CODE_FASTMARCHCDTW_H
