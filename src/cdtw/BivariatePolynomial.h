#ifndef TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H

#include <array>
#include <iostream>
#include "Polynomial.h"

template<size_t D>
struct BivariatePolynomial
{
    std::array<std::array<double, D + 1>, D + 1> coefficients;

    BivariatePolynomial() = default;
    explicit BivariatePolynomial(const std::array<std::array<double, D + 1>, D + 1>& coefficients) :
        coefficients(coefficients) {}

    friend std::ostream& operator<<(std::ostream& os, const BivariatePolynomial& p) {
        for (int dx = D; dx >= 0; --dx) {
            for (int dy = D; dy >= 0; --dy) {
                if (dx + dy > D) continue;

                if (dx < D) os << " + ";
                os << p.coefficients[dx][dy];
                if (dx > 0) os << "*x";
                if (dx > 1) os << "^" << dx;
                if (dy > 0) os << "*y";
                if (dy > 1) os << "^" << dy;
            }
        }
        return os;
    }

    BivariatePolynomial<D - 1> partial_derivative_x() const {
        BivariatePolynomial<D - 1> result;
        for (size_t dx = 0; dx < D; ++dx) {
            for (size_t dy = 0; dy < D; ++dy) {
                result.coefficients[dx][dy] = coefficients[dx + 1][dy] * (dx + 1);
            }
        }
        return result;
    }

    BivariatePolynomial<D - 1> partial_derivative_y() const {
        BivariatePolynomial<D - 1> result;
        for (size_t dx = 0; dx < D; ++dx) {
            for (size_t dy = 0; dy < D; ++dy) {
                result.coefficients[dx][dy] = coefficients[dx][dy + 1] * (dy + 1);
            }
        }
        return result;
    }

    /**
     * Sum this function with a function in x.
     * Returns f(x, y) + g(x).
     */
    BivariatePolynomial<D> add_x(const Polynomial<D>& g) const {
        auto result = *this;
        for (size_t dx = 0; dx <= D; ++dx) {
            result.coefficients[dx][0] += g.coefficients[dx];
        }
        return result;
    }

    BivariatePolynomial<D> translate_xy(double x, double y) const;

    //
    // Arithmetic operators
    //

    // BivariatePolynomial + BivariatePolynomial (of same degree)
    BivariatePolynomial<D>& operator+=(const BivariatePolynomial<D>& rhs) {
        for (size_t dx = 0; dx <= D; ++dx) {
            for (size_t dy = 0; dy <= D; ++dy) {
                coefficients[dx][dy] += rhs.coefficients[dx][dy];
            }
        }
        return *this;
    }
    friend BivariatePolynomial<D> operator+(BivariatePolynomial<D> lhs, const BivariatePolynomial<D>& rhs) {
        lhs += rhs;
        return lhs;
    }
};

/**
 * Find the functions g(y) where f(g(y), y) = 0.
 *
 * @param f
 * @return
 */
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<1>& f) {
    auto c00 = f.coefficients[0][0];
    auto c01 = f.coefficients[0][1];
    auto c10 = f.coefficients[1][0];

    if (c10 == 0) {
        return {};
    }

    return {
        Polynomial<1>({-c00 / c10, -c01 / c10})
    };
}

template<>
BivariatePolynomial<2> BivariatePolynomial<2>::translate_xy(double cx, double cy) const {
    auto c00 = coefficients[0][0];
    auto c01 = coefficients[0][1];
    auto c02 = coefficients[0][2];
    auto c10 = coefficients[1][0];
    auto c20 = coefficients[2][0];
    auto c11 = coefficients[1][1];

    return BivariatePolynomial<2>({{
        {{c00 - c10*cx + c20*cx*cx - c01*cy + c11*cx*cy + c02*cy*cy, c01 - c11*cx - 2*c02*cy, c02}},
        {{c10 - 2*c20*cx - c11*cy, c11, 0}},
        {{c20, 0, 0}},
    }});
}

template<size_t D>
std::ostream& operator<<(std::ostream& os, const BivariatePolynomial<D>& p) {
    for (size_t dx = 0; dx < D; ++dx) {
        for (size_t dy = 0; dy < D; ++dy) {
            os << p.coefficients[dx][dy] << "*x^" << dx << "*y^" << dy;
        }
    }
    return os;
}

#endif //TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
