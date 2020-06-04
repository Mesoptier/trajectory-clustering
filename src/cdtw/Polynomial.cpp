#include "Polynomial.h"

#include <unsupported/Eigen/Polynomials>

//
// find_roots
//

template<>
std::vector<double> find_roots<1>(const Polynomial<1>& f) {
    double c0 = f.coefficients[0];
    double c1 = f.coefficients[1];

    if (c1 == 0) {
        return {};
    }
    return {-c0 / c1};
}

template<>
std::vector<double> find_roots<2>(const Polynomial<2>& f) {
    double c0 = f.coefficients[0];
    double c1 = f.coefficients[1];
    double c2 = f.coefficients[2];

    if (c2 == 0) {
        return find_roots(Polynomial<1>({c0, c1}));
    }

    double discriminant = c1 * c1 - 4 * c2 * c0;
    if (approx_zero(discriminant)) {
        return {-c1 / (2 * c2)};
    }
    if (discriminant < 0) {
        // Complex roots
        return {};
    }

    return {
        (-c1 + sqrt(discriminant)) / (2 * c2),
        (-c1 - sqrt(discriminant)) / (2 * c2),
    };
}

template<>
std::vector<double> find_roots<3>(const Polynomial<3>& f) {
    double d = f.coefficients[0];
    double c = f.coefficients[1];
    double b = f.coefficients[2];
    double a = f.coefficients[3];

    if (approx_zero(a)) {
        return find_roots(Polynomial<2>({d, c, b}));
    }

    Eigen::PolynomialSolver<double, 3> solver;
    Eigen::Array4d poly(d, c, b, a);
    solver.compute(poly);

    std::vector<double> roots;
    solver.realRoots(roots);
    return roots;
}

//
// translate_xy
//

template<>
Polynomial<1> Polynomial<1>::translate_xy(double cx, double cy) const {
    auto c0 = coefficients[0];
    auto c1 = coefficients[1];

    return Polynomial<1>({c0 + cy - c1 * cx, c1});
}

template<>
Polynomial<2> Polynomial<2>::translate_xy(double cx, double cy) const {
    auto c0 = coefficients[0];
    auto c1 = coefficients[1];
    auto c2 = coefficients[2];

    return Polynomial<2>({c0 - c1 * cx + c2 * cx * cx + cy, c1 - 2 * c2 * cx, c2});
}