

#include <array>
#include <iostream>
#include "BivariatePolynomial.h"

// template<>
// std::ostream& operator<<(std::ostream& os, const BivariatePolynomial<1>& p) {
//         for (int dx = 1; dx >= 0; --dx) {
//             for (int dy = 1; dy >= 0; --dy) {
//                 if (dx + dy > 1) continue;

//                 if (dx < 1) os << " + ";
//                 os << p.coefficients[dx][dy];
//                 if (dx > 0) os << "*x";
//                 if (dx > 1) os << "^" << dx;
//                 if (dy > 0) os << "*y";
//                 if (dy > 1) os << "^" << dy;
//             }
//         }
//         return os;
//     }

// template<>
// std::ostream& operator<<(std::ostream& os, const BivariatePolynomial<2>& p) {
//         for (int dx = 2; dx >= 0; --dx) {
//             for (int dy = 2; dy >= 0; --dy) {
//                 if (dx + dy > 2) continue;

//                 if (dx < 2) os << " + ";
//                 os << p.coefficients[dx][dy];
//                 if (dx > 0) os << "*x";
//                 if (dx > 1) os << "^" << dx;
//                 if (dy > 0) os << "*y";
//                 if (dy > 1) os << "^" << dy;
//             }
//         }
//         return os;
//     }

    template<>
    double BivariatePolynomial<2>::operator()(Point p) {
        
        double result = 0;

        for (int i = 0; i <= 2; ++i)
            for (int j = 0; j <= 2; ++j) {
                result += coefficients[i][j] * pow(p.x, i) * pow(p.y, j);
            }

        return result;
    }

    template<>
    double BivariatePolynomial<1>::operator()(Point p) {
        
        double result = 0;

        for (int i = 0; i <= 1; ++i)
            for (int j = 0; j <= 1; ++j) {
                result += coefficients[i][j] * pow(p.x, i) * pow(p.y, j);
            }

        return result;
    }

        template<>
        BivariatePolynomial<1> BivariatePolynomial<2>::partial_derivative_x() const {
        BivariatePolynomial<1> result;
        for (size_t dx = 0; dx < 2; ++dx) {
            for (size_t dy = 0; dy < 2; ++dy) {
                result.coefficients[dx][dy] = coefficients[dx + 1][dy] * (dx + 1);
            }
        }
        return result;
    }

    template<>
    BivariatePolynomial<1> BivariatePolynomial<2>::partial_derivative_y() const {
    BivariatePolynomial<1> result;
        for (size_t dx = 0; dx < 2; ++dx) {
            for (size_t dy = 0; dy < 2; ++dy) {
                result.coefficients[dx][dy] = coefficients[dx][dy + 1] * (dy + 1);
            }
        }
        return result;
    }


    template<>
    BivariatePolynomial<1> BivariatePolynomial<1>::add_x(const Polynomial<1>& g) const {
        auto result = *this;
        for (size_t dx = 0; dx <= 1; ++dx) {
            result.coefficients[dx][0] += g.coefficients[dx];
        }
        return result;
    }

    template<>
    BivariatePolynomial<2> BivariatePolynomial<2>::add_x(const Polynomial<2>& g) const {
        auto result = *this;
        for (size_t dx = 0; dx <= 2; ++dx) {
            result.coefficients[dx][0] += g.coefficients[dx];
        }
        return result;
    }

    

    template<>
    Polynomial<1> BivariatePolynomial<1>::slice_at_y(double c) const {
        Polynomial<1> result;
        for (size_t dx = 0; dx <= 1; ++dx) {
            for (size_t dy = 0; dy <= 1; ++dy) {
                if (dy == 0) {
                    result.coefficients[dx] += coefficients[dx][dy];
                } else {
                    result.coefficients[dx] += coefficients[dx][dy] * std::pow(c, dy);
                }
            }
        }
        return result;
    }

    template<>
    Polynomial<2> BivariatePolynomial<2>::slice_at_y(double c) const {
        Polynomial<2> result;
        for (size_t dx = 0; dx <= 2; ++dx) {
            for (size_t dy = 0; dy <= 2; ++dy) {
                if (dy == 0) {
                    result.coefficients[dx] += coefficients[dx][dy];
                } else {
                    result.coefficients[dx] += coefficients[dx][dy] * std::pow(c, dy);
                }
            }
        }
        return result;
    }

    //
    // Arithmetic operators
    //

    // BivariatePolynomial + BivariatePolynomial (of same degree)
    template<>
    BivariatePolynomial<1>& BivariatePolynomial<1>::operator+=(const BivariatePolynomial<1>& rhs) {
        for (size_t dx = 0; dx <= 1; ++dx) {
            for (size_t dy = 0; dy <= 1; ++dy) {
                coefficients[dx][dy] += rhs.coefficients[dx][dy];
            }
        }
        return *this;
    }

    template<>
    BivariatePolynomial<2>& BivariatePolynomial<2>::operator+=(const BivariatePolynomial<2>& rhs) {
        for (size_t dx = 0; dx <= 2; ++dx) {
            for (size_t dy = 0; dy <= 2; ++dy) {
                coefficients[dx][dy] += rhs.coefficients[dx][dy];
            }
        }
        return *this;
    }

    template<>
    BivariatePolynomial<1>& BivariatePolynomial<1>::operator*=(double a) {
        for (size_t dx = 0; dx <= 1; ++dx) {
            for (size_t dy = 0; dy <= 1; ++dy) {
                coefficients[dx][dy] *= a;
            }
        }
        return *this;
    }
    template<>
    BivariatePolynomial<2>& BivariatePolynomial<2>::operator*=(double a) {
        for (size_t dx = 0; dx <= 2; ++dx) {
            for (size_t dy = 0; dy <= 2; ++dy) {
                coefficients[dx][dy] *= a;
            }
        }
        return *this;
    } 

    template<>
    BivariatePolynomial<1> BivariatePolynomial<1>::operator+(const BivariatePolynomial<1>& rhs) {
        BivariatePolynomial<1> aux = *this;
        aux += rhs;
        return aux;
    }

    template<>
    BivariatePolynomial<2> BivariatePolynomial<2>::operator+(const BivariatePolynomial<2>& rhs) {
        BivariatePolynomial<2> aux = *this;
        aux += rhs;
        return aux;
    }


// template<size_t D>
// struct BivariatePolynomial
// {
//     std::array<std::array<double, D + 1>, D + 1> coefficients;

//     BivariatePolynomial() = default;
//     explicit BivariatePolynomial(const std::array<std::array<double, D + 1>, D + 1>& coefficients) :
//         coefficients(coefficients) {}


// };

/**
 * Find the functions g(y) where f(g(y), y) = 0.
 *
 * @param f
 * @return
 */
template<size_t D>
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<D>& f);

template<>
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
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<2>& f) {
    auto c00 = f.coefficients[0][0];
    auto c01 = f.coefficients[0][1];
    auto c02 = f.coefficients[0][2];
    auto c10 = f.coefficients[1][0];
    auto c11 = f.coefficients[1][1];
    auto c20 = f.coefficients[2][0];

    if (approx_zero(c11) && approx_zero(c01) && approx_zero(c02)) {
        const std::vector<double> roots = find_roots(Polynomial<2>({c00, c10, c20}));
        std::vector<Polynomial<1>> result;
        for (const auto root : roots) {
            result.push_back(Polynomial<1>({root, 0}));
        }
        return result;
    }

    if (approx_zero(c00) && approx_zero(c10) && approx_zero(c20)) {
        return {
            Polynomial<1>({-c01 / c11, -c02 / c11}),
        };
    }

    if (c11 == 4 && c20 == -2 && c02 == -2 && c10 == -c01) {
        double d = 8 * c00 + c10 * c10;
        if (approx_zero(d)) {
            return {
                Polynomial<1>({(c10 - std::sqrt(d)) / 4, 1}),
            };
        }
        if (d < 0) {
            return {};
        }
        return {
            Polynomial<1>({(c10 - std::sqrt(d)) / 4, 1}),
            Polynomial<1>({(c10 + std::sqrt(d)) / 4, 1}),
        };
    }

    // TODO: Maybe use sampling to find a few roots and return a set of linear functions through those points?

    return {};
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

// TODO: Generalize embed_x

template<>
Polynomial<2> BivariatePolynomial<2>::embed_x(const Polynomial<1>& g) const {
    auto c00 = coefficients[0][0];
    auto c01 = coefficients[0][1];
    auto c02 = coefficients[0][2];
    auto c10 = coefficients[1][0];
    auto c20 = coefficients[2][0];
    auto c11 = coefficients[1][1];

    auto g0 = g.coefficients[0];
    auto g1 = g.coefficients[1];

    return Polynomial<2>({
        c00 + c10*g0 + c20*g0*g0,
        c01 + c11*g0 + c10*g1 + 2*c20*g0*g1,
        c02 + c11*g1 + c20*g1*g1,
    });
}

template<>
Polynomial<3> BivariatePolynomial<3>::embed_x(const Polynomial<1>& g) const {
    auto c00 = coefficients[0][0];
    auto c01 = coefficients[0][1];
    auto c02 = coefficients[0][2];
    auto c03 = coefficients[0][3];
    auto c10 = coefficients[1][0];
    auto c11 = coefficients[1][1];
    auto c12 = coefficients[1][2];
    auto c20 = coefficients[2][0];
    auto c21 = coefficients[2][1];
    auto c30 = coefficients[3][0];

    auto g0 = g.coefficients[0];
    auto g1 = g.coefficients[1];

    return Polynomial<3>({
        c00 + c10*g0 + c20*g0*g0 + c30*g0*g0*g0,
        c01 + c11*g0 + c21*g0*g0 + c10*g1 + 2*c20*g0*g1 + 3*c30*g0*g0*g1,
        c02 + c12*g0 + c11*g1 + 2*c21*g0*g1 + c20*g1*g1 + 3*c30*g0*g1*g1,
        c03 + c12*g1 + c21*g1*g1 + c30*g1*g1*g1,
    });
}


