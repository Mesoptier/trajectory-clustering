#ifndef TRAJECTORY_CLUSTERING_POLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_POLYNOMIAL_H

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include "../util.h"

template<size_t D>
struct Polynomial
{
    std::array<double, D + 1> coefficients;

    Polynomial() : coefficients() {}
    explicit Polynomial(const std::array<double, D + 1>& coefficients) : coefficients(coefficients) {}

    /**
     * Evaluate the polynomial at x.
     * @param x
     * @return
     */
    double operator()(double x) const {
        double y = 0;
        for (size_t d = 0; d <= D; ++d) {
            y += coefficients[d] * std::pow(x, d);
        }
        return y;
    }

    /**
     * Get the first order derivative.
     * @return
     */
    Polynomial<D> derivative() const {
        Polynomial<D> result;
        for (size_t d = 0; d < D; ++d) {
            result.coefficients[d] = coefficients[d + 1] * (d + 1);
        }
        return result;
    }

    bool changes_sign_at(double x) const {
        if (!approx_zero((*this)(x))) {
            return false;
        }
        return !this->derivative().changes_sign_at(x);
    }

    Polynomial<D> translate_xy(double cx, double cy) const;

    //
    // Arithmetic operators
    //

    // Polynomial + Polynomial (of same degree)
    Polynomial<D>& operator+=(const Polynomial<D>& rhs) {
        for (size_t d = 0; d <= D; ++d) {
            coefficients[d] += rhs.coefficients[d];
        }
        return *this;
    }
    friend Polynomial<D> operator+(Polynomial<D> lhs, const Polynomial<D>& rhs) {
        lhs += rhs;
        return lhs;
    }

    // Polynomial + constant
    Polynomial<D>& operator+=(double c) {
        coefficients[0] += c;
        return *this;
    }
    friend Polynomial<D> operator+(Polynomial<D> f, double c) {
        f += c;
        return f;
    }

    Polynomial<D>& operator-=(const Polynomial<D>& rhs) {
        for (size_t d = 0; d <= D; ++d) {
            coefficients[d] -= rhs.coefficients[d];
        }
        return *this;
    }
    friend Polynomial<D> operator-(Polynomial<D> lhs, const Polynomial<D>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    //
    // Equality and relational operators
    //

    bool operator==(const Polynomial<D>& rhs) const {
        return coefficients == rhs.coefficients;
    }

    bool operator!=(const Polynomial<D>& rhs) const {
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
    struct CompareAt
    {
        double x;
        explicit CompareAt(double x) : x(x) {}
        bool operator()(Polynomial<D> f, Polynomial<D> g) const {
            size_t d = D;
            while (d > 0 && approx_equal(f(x), g(x))) {
                f = f.derivative();
                g = g.derivative();
                --d;
            }
            return f(x) < g(x);
        }
    };
};

template<size_t D>
std::ostream& operator<<(std::ostream& os, const Polynomial<D>& p) {
    for (int d = D; d >= 0; --d) {
        if (d < D) os << " + ";
        os << p.coefficients[d];
        if (d > 0) os << "*x";
        if (d > 1)os << "^" << d;
    }
    return os;
}

/**
 * Find the values for x where f(x) = 0.
 *
 * @param f
 * @return
 */
template<size_t D>
std::vector<double> find_roots(const Polynomial<D>& f);

std::vector<double> find_roots(const Polynomial<1>& f) {
    double c0 = f.coefficients[0];
    double c1 = f.coefficients[1];

    if (c1 == 0) {
        return {};
    }
    return {-c0 / c1};
}

/**
 * Find the values for x where f(x) = 0.
 *
 * @param f
 * @return
 */
std::vector<double> find_roots(const Polynomial<2>& f) {
    double c0 = f.coefficients[0];
    double c1 = f.coefficients[1];
    double c2 = f.coefficients[2];

    if (c2 == 0) {
        if (c1 == 0) {
            return {};
        }
        return {-c0 / c1};
    }

    double det = c1 * c1 - 4 * c2 * c0;
    if (det == 0) {
        return {-c1 / (2 * c2)};
    }
    if (det < 0) {
        throw std::runtime_error("Complex roots");
    }

    return {
        (-c1 + sqrt(det)) / (2 * c2),
        (-c1 - sqrt(det)) / (2 * c2),
    };
}

/**
 * Find the values for x where f(x) = g(x).
 *
 * @tparam D
 * @param f
 * @param g
 * @return
 */
template<size_t D>
std::vector<double> find_intersections(const Polynomial<D>& f, const Polynomial<D>& g) {
    return find_roots(f - g);
}

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

#endif //TRAJECTORY_CLUSTERING_POLYNOMIAL_H
