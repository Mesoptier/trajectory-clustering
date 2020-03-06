#ifndef TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H

#include <array>
#include <iostream>

template<size_t D>
struct BivariatePolynomial {
    std::array<std::array<double, D + 1>, D + 1> coefficients;

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
};

#endif //TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
