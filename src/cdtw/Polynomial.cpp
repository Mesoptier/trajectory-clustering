#include "Polynomial.h"

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

    // Depressed cubic equation
    // t^3 + pt + q = 0
    double shift = b / (3*a); // t = x + shift
    double p = (3*a*c - b*b) / (3*a*a);
    double q = (2*b*b*b - 9*a*b*c + 27*a*a*d) / (27*a*a*a);

    double D = -(4*p*p*p + 27*q*q);

    if (approx_zero(D)) {
        if (approx_zero(p)) { // 1 triple root
            double t1 = 0;
            return { t1 - shift };
        }

        // 1 simple root + 1 double root
        double t1 = (3*q) / p;
        double t2 = (-3*q) / (2*p);
        return { t1 - shift, t2 - shift };
    }

    if (D < 0) { // 1 real root + 2 complex roots
        double term1 = std::sqrt((q*q)/4 + (p*p*p)/27);
        double t1 = std::cbrt(-q/2 + term1) + std::cbrt(-q/2 - term1);
        return { t1 - shift };
    }

    // 3 real roots
    double term1 = 2 * std::sqrt(-p/3);
    double term2 = 1./3 * std::acos((3*q) / (2*p) * std::sqrt(-3/p));
    std::vector<double> roots(3);

    for (size_t k = 0; k < 3; ++k) {
        roots[k] = term1 * std::cos(term2 - (2 * M_PIl * k) / 3) - shift;
    }

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