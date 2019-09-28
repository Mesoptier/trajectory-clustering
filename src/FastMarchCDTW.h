#ifndef CODE_FASTMARCHCDTW_H
#define CODE_FASTMARCHCDTW_H

#include <armadillo>
#include <utility>
#include <vector>
#include "Curve.h"

typedef std::pair<int, int> Point;

inline Point operator +(const Point& lhs, const Point& rhs) {
    return {lhs.first + rhs.first, lhs.second + rhs.second};
}

class FastMarchCDTW {

    // User options
    double h;
    int imageNorm = 2;
    int paramNorm = 1;

    // Own variables (initialized in ::compute)
    arma::mat f_mat;
    arma::mat u_mat;
    unsigned int n_rows = 0;
    unsigned int n_cols = 0;
    double hi = 0;
    double hj = 0;


public:
    FastMarchCDTW(double h, int imageNorm, int paramNorm);

    double compute(const Curve<double>& curve1, const Curve<double>& curve2, bool saveMatrices = false);


private:
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

    bool inBounds(Point point, unsigned int n_rows, unsigned int n_cols);

    double eikonalUpdate(int i, int j);
};


#endif //CODE_FASTMARCHCDTW_H
