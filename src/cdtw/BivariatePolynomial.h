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

    return {
        Polynomial<1>({-c00 / c10, -c01 / c10})
    };
}

#endif //TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
