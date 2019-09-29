#ifndef CODE_FASTMARCHINTEGRALFRECHET_H
#define CODE_FASTMARCHINTEGRALFRECHET_H

#include <armadillo>
#include <utility>
#include <vector>
#include "Curve.h"

typedef std::pair<int, int> Point;

inline Point operator +(const Point& lhs, const Point& rhs) {
    return {lhs.first + rhs.first, lhs.second + rhs.second};
}

class FastMarchIntegralFrechet {

    // User options
    Curve<double> curve1;
    Curve<double> curve2;
    int imageNorm;
    int paramNorm;

    // Own variables
    unsigned int n_rows = 0;
    unsigned int n_cols = 0;
    double hi = 0;
    double hj = 0;
    arma::mat f_mat;
    arma::mat u_mat;

    arma::mat matching;
    arma::mat center;

public:
    FastMarchIntegralFrechet(const Curve<double>& curve1, const Curve<double>& curve2, double h, int imageNorm = 2, int paramNorm = 1);

    double computeDistance();

    void computeMatching(double stepSize = 0.01, int maxIterations = 10000);

    void computeCenter(double ratio);

    void save();

    const arma::mat& getCenter() const;

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

    bool inBounds(Point point);

    double eikonalUpdate(int i, int j);
};


#endif //CODE_FASTMARCHINTEGRALFRECHET_H
