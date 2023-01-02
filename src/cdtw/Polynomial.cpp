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

 // Polynomial + Polynomial (of same degree)
    // template<size_t D>
    // Polynomial<D>& Polynomial<D>::operator+=(const Polynomial<D>& rhs) {
    //     for (size_t d = 0; d <= D; ++d) {
    //         coefficients[d] += rhs.coefficients[d];
    //     }
    //     return *this;
    // }

    template<>
    Polynomial<1> Polynomial<1>::operator+(const Polynomial<1>& rhs) const {
        Polynomial<1> aux = *this;
        aux += rhs;
        return aux;
    }

    template<>
    Polynomial<2> Polynomial<2>::operator+(const Polynomial<2>& rhs) const {
        Polynomial<2> aux = *this;
        aux += rhs;
        return aux;
    }

    // Polynomial + constant
    // template<size_t D>
    // Polynomial<D>& Polynomial<D>::operator+=(double c) {
    //     coefficients[0] += c;
    //     return *this;
    // }


    template<>
    Polynomial<1> Polynomial<1>::operator+(double c) {
        *this += c;
        return *this;
    }

    template<>
    Polynomial<2> Polynomial<2>::operator+(double c) {
        *this += c;
        return *this;
    }


    template<>
    Polynomial<1>& Polynomial<1>::operator-=(const Polynomial<1>& rhs) {
        for (size_t d = 0; d <= 1; ++d) {
            coefficients[d] -= rhs.coefficients[d];
        }
        return *this;
    }

    template<>
    Polynomial<2>& Polynomial<2>::operator-=(const Polynomial<2>& rhs) {
        for (size_t d = 0; d <= 2; ++d) {
            coefficients[d] -= rhs.coefficients[d];
        }
        return *this;
    }

    template<>
    Polynomial<1> Polynomial<1>::operator-(const Polynomial<1>& rhs) const {
        Polynomial<1> aux = *this;
        aux -= rhs;
        return aux;
    }

    template<>
    Polynomial<2> Polynomial<2>::operator-(const Polynomial<2>& rhs) const {
        Polynomial<2> aux = *this;
        aux -= rhs;
        return aux;
    }

    //
    // Equality and relational operators
    //

    template<>
    bool Polynomial<1>::operator==(const Polynomial<1>& rhs) const {
        return coefficients == rhs.coefficients;
    }


    template<>
    bool Polynomial<2>::operator==(const Polynomial<2>& rhs) const {
        for (int i = 0; i < coefficients.size(); ++i)
            if (!approx_equal(coefficients[i], rhs.coefficients[i]))
                return false;
        return true;
        return coefficients == rhs.coefficients;
    }

    template<>
    bool Polynomial<1>::operator!=(const Polynomial<1>& rhs) const {
        return !(rhs == *this);
    }

    template<>
    bool Polynomial<2>::operator!=(const Polynomial<2>& rhs) const {
        return !(rhs == *this);
    }

    /**
     * Compare two polynomials of the same degree lexicographically first by value, then by value of first derivative,
     * then second derivative, etc. at a given x.
     *
     * That is, f(x) is considered smaller than g(x) if:
     *      f(x) < g(x)
     *   || (f(x) = g(x) && f'(x) < g'(x))
     *   || (f(x) = g(x) && f'(x) = g'(x) && f''(x) < g''(x))
     *   || ...
     */

    template<>
    std::string Polynomial<1>::to_string() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    template<>
    std::string Polynomial<2>::to_string() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }