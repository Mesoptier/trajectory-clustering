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
    if (approx_zero(f.coefficients.back())) {
        return find_roots(change_degree<1>(f));
    }

    Eigen::PolynomialSolver<double, 2> solver;
    solver.compute(Eigen::Array3d(f.coefficients.data()));

    std::vector<double> roots;
    solver.realRoots(roots);
    return roots;
}

template<>
std::vector<double> find_roots<3>(const Polynomial<3>& f) {
    if (approx_zero(f.coefficients.back())) {
        return find_roots(change_degree<2>(f));
    }

    Eigen::PolynomialSolver<double, 3> solver;
    solver.compute(Eigen::Array4d(f.coefficients.data()));

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