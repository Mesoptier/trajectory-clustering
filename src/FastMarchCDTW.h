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

    double compute(
        const Curve<double>& curve1,
        const Curve<double>& curve2,
        double h,
        int imageNorm = 2,
        int paramNorm = 2,
        bool saveMatrices = false
    );

    bool inBounds(Point point, unsigned int n_rows, unsigned int n_cols);

    double eikonalUpdate(const arma::mat& u_mat, const arma::mat& f_mat, int i, int j, double hi, double hj, int norm);
}


#endif //CODE_FASTMARCHCDTW_H
