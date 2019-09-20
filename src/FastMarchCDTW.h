#ifndef CODE_FASTMARCHCDTW_H
#define CODE_FASTMARCHCDTW_H

#include <armadillo>
#include <utility>
#include <vector>
#include "Curve.h"

namespace FastMarchCDTW {

    typedef std::pair<int, int> Point;

    inline Point operator +(const Point& lhs, const Point& rhs) {
        return {lhs.first + rhs.first, lhs.second + rhs.second};
    }

    struct Node {
        Node(Point point, double cost) : point(point), cost(cost) {}

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

    double compute(const Curve<double>& curve1, const Curve<double>& curve2, double h, bool saveMatrices = false);

    bool inBounds(Point point, unsigned int n_rows, unsigned int n_cols);
}


#endif //CODE_FASTMARCHCDTW_H
