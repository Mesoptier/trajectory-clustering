#ifndef TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H

#include <ostream>
#include "PiecewisePolynomial.h"
#include "BivariatePolynomial.h"

// TODO: Rename to BivariatePolynomialPiece

template<size_t D>
struct ConstrainedBivariatePolynomial
{
    // f(x, y)
    BivariatePolynomial<D> f;
    // y is valid in this interval
    Interval y_interval;
    // x is valid when g(y) <= x <= h(y) for each g in left_constraints and h in right_constraints
    std::vector<Polynomial<1>> left_constraints;
    std::vector<Polynomial<1>> right_constraints;

    /**
     * Sum this function with a function in x.
     * Returns f(x, y) + g(x) with updated constraints.
     */
    ConstrainedBivariatePolynomial<D> add_x(const PolynomialPiece<D>& g) const;

    ConstrainedBivariatePolynomial<D> translate_xy(double cx, double cy) const {
        auto result = *this;
        result.f = result.f.translate_xy(cx, cy);
        result.y_interval.min += cy;
        result.y_interval.max += cy;
        for (auto& c : result.left_constraints) {
            c = c.translate_xy(cy, cx);
        }
        for (auto& c : result.right_constraints) {
            c = c.translate_xy(cy, cx);
        }
        return result;
    }

    PolynomialPiece<D> slice_at_y(double c) const {
        return {interval_at_y(c), f.slice_at_y(c)};
    }

    Interval interval_at_y(double y) const {
        double min = -std::numeric_limits<double>::infinity();
        double max = std::numeric_limits<double>::infinity();

        for (const auto& lc : left_constraints) {
            min = std::max(min, lc(y));
        }
        for (const auto& rc : right_constraints) {
            max = std::min(max, rc(y));
        }

        return { min, max };
    }

    friend std::ostream& operator<<(std::ostream& os, const ConstrainedBivariatePolynomial& polynomial) {
        os << std::fixed;
        os << "{(" << polynomial.f << "), ";
        os << "((" << polynomial.y_interval << ") /. x->y)";
        for (auto& c : polynomial.left_constraints) {
            os << " && ((" << c << ") /. x->y) <= x";
        }
        for (auto& c : polynomial.right_constraints) {
            os << " && ((" << c << ") /. x->y) >= x";
        }
        os << "}";
        return os;
    }
};

template<size_t D>
ConstrainedBivariatePolynomial<D> ConstrainedBivariatePolynomial<D>::add_x(const PolynomialPiece<D>& g) const {
    auto result = *this;
    result.f = result.f.add_x(g.polynomial);
    result.left_constraints.push_back(Polynomial<1>({g.interval.min, 0}));
    result.right_constraints.push_back(Polynomial<1>({g.interval.max, 0}));
    return result;
}

#endif //TRAJECTORY_CLUSTERING_CONSTRAINEDBIVARIATEPOLYNOMIAL_H
