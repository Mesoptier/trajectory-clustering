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

    double operator()(Point p);

    BivariatePolynomial<D - 1> partial_derivative_x() const;

    BivariatePolynomial<D - 1> partial_derivative_y() const;

    /**
     * Sum this function with a function in x.
     * Returns f(x, y) + g(x).
     */
    BivariatePolynomial<D> add_x(const Polynomial<D>& g) const;

    /**
     * Returns f(g(y), y)
     */
    Polynomial<D> embed_x(const Polynomial<1>& g) const;

    /**
     * Returns g(x) = f(x, c)
     */
    Polynomial<D> slice_at_y(double c) const;

    BivariatePolynomial<D> translate_xy(double x, double y) const;

    //
    // Arithmetic operators
    //

    // BivariatePolynomial + BivariatePolynomial (of same degree)
    BivariatePolynomial<D>& operator+=(const BivariatePolynomial<D>& rhs);

    BivariatePolynomial<D>& operator*=(double a);

    BivariatePolynomial<D> operator+(const BivariatePolynomial<D>& rhs);

    BivariatePolynomial<D> operator*(double a);

    // friend std::ostream& operator<<(std::ostream& os, const BivariatePolynomial<D>& p);


};


/**
 * Find the functions g(y) where f(g(y), y) = 0.
 *
 * @param f
 * @return
 */
template<size_t D>
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<D>& f);

template<>
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<1>& f);

template<>
std::vector<Polynomial<1>> find_roots_y(const BivariatePolynomial<2>& f);

template<>
BivariatePolynomial<2> BivariatePolynomial<2>::translate_xy(double cx, double cy) const;

// TODO: Generalize embed_x

template<>
Polynomial<2> BivariatePolynomial<2>::embed_x(const Polynomial<1>& g) const;

template<>
Polynomial<3> BivariatePolynomial<3>::embed_x(const Polynomial<1>& g) const;

template<size_t D>
std::ostream& operator<<(std::ostream& os, const BivariatePolynomial<D>& p);

#endif //TRAJECTORY_CLUSTERING_BIVARIATEPOLYNOMIAL_H
