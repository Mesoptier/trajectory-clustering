#ifndef TRAJECTORY_CLUSTERING_POLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_POLYNOMIAL_H

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
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
            if (d == 0) {
                y += coefficients[d];
            } else if (d == 1) {
                y += coefficients[d] * x;
            } else {
                y += coefficients[d] * std::pow(x, d);
            }
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

    /**
     * Test whether the function changes sign at x.
     * @param x
     * @return
     */
    bool changes_sign_at(double x) const {
        // For a function to change sign at x:
        //  1. It must equal zero at x
        //  2. Its derivative must nog change sign at x

        Polynomial<D> f = *this;
        bool state = false;

        for (size_t d = D; d > 0; --d) {
            if (!approx_zero(f(x))) {
                return state;
            }
            state = !state;
            f = f.derivative();
        }

        return state;
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
            if (approx_equal(x, 0.000170711)) {
                std::cout << f(x) << std::endl;
                std::cout << g(x) << std::endl;
                std::cout << x << std::endl;
            }
            return f(x) < g(x);
        }
    };

    inline std::string to_string() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }
};

template<size_t D>
bool approx_equal(const Polynomial<D>& a, const Polynomial<D>& b, double tol = ABS_TOL) {
    return approx_equal(a.coefficients, b.coefficients, tol);
}

template<size_t NewD, size_t D>
Polynomial<NewD> change_degree(const Polynomial<D>& f) {
    Polynomial<NewD> result;
    std::copy(f.coefficients.begin(), f.coefficients.begin() + std::min(D, NewD) + 1, result.coefficients.begin());
    return result;
}

template<size_t D>
std::ostream& operator<<(std::ostream& os, const Polynomial<D>& p) {
    os << std::fixed;
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

#endif //TRAJECTORY_CLUSTERING_POLYNOMIAL_H
