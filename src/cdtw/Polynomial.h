#ifndef TRAJECTORY_CLUSTERING_POLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_POLYNOMIAL_H

#include <array>
#include <vector>
#include <cmath>
#include <iostream>

template<size_t D>
struct Polynomial
{
    std::array<double, D + 1> coefficients;

    Polynomial() = default;
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

    //
    // Arithmetic operators
    //

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
            while (d > 0 && f(x) == g(x)) {
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
std::vector<double> find_roots(const Polynomial<1>& f) {
    double a = f.coefficients[1];
    double b = f.coefficients[0];

    if (a == 0) {
        if (b == 0) {
            throw std::runtime_error("f is always 0");
        }
        return {};
    }
    return {-b / a};
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


#endif //TRAJECTORY_CLUSTERING_POLYNOMIAL_H
