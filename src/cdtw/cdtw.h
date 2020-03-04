#pragma once

#include <vector>
#include <array>
#include <ostream>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>

struct Interval {
    double min;
    double max;

    bool contains(double x) const {
        return min <= x && x <= max;
    }

    friend std::ostream& operator<<(std::ostream& os, const Interval& interval) {
        os << "[" << interval.min << ", " << interval.max << "]";
        return os;
    }
};

template<size_t D>
struct Polynomial {
    std::array<double, D + 1> coefficients;

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
        for (int d = D; d >= 0; --d) {
            if (d < D) os << " + ";
            os << p.coefficients[d];
            if (d > 0) os << "*x";
            if (d > 1)os << "^" << d;
        }
        return os;
    }

    double operator()(double x) const {
        double y = 0;
        for (size_t d = 0; d <= D; ++d) {
            y += coefficients[d] * std::pow(x, d);
        }
        return y;
    }

    Polynomial derivative() const {
        Polynomial result;
        for (size_t d = 0; d < D; ++d) {
            result.coefficients[d] = coefficients[d + 1] * (d + 1);
        }
        return result;
    }

    //
    // Arithmetic operators
    //

    Polynomial& operator+=(const Polynomial& rhs) {
        for (size_t d = 0; d <= D; ++d) {
            coefficients[d] += rhs.coefficients[d];
        }
        return *this;
    }
    friend Polynomial operator+(Polynomial lhs, const Polynomial& rhs) {
        lhs += rhs;
        return lhs;
    }

    Polynomial& operator-=(const Polynomial& rhs) {
        for (size_t d = 0; d <= D; ++d) {
            coefficients[d] -= rhs.coefficients[d];
        }
        return *this;
    }
    friend Polynomial operator-(Polynomial lhs, const Polynomial& rhs) {
        lhs -= rhs;
        return lhs;
    }

    //
    // Equality and relational operators
    //

    bool operator==(const Polynomial& rhs) const {
        return coefficients == rhs.coefficients;
    }

    bool operator!=(const Polynomial& rhs) const {
        return !(rhs == *this);
    }

    struct CompareAt {
        double x;
        explicit CompareAt(double x) : x(x) {}
        bool operator()(Polynomial f, Polynomial g) const {
            size_t d = D;
            while (f(x) == g(x) && d > 0) {
                f = f.derivative();
                g = g.derivative();
                --d;
            }
            return f(x) < g(x);
        }
    };
};

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

template<size_t D>
std::vector<double> find_intersections(const Polynomial<D>& f, const Polynomial<D>& g) {
    return find_roots(f - g);
}

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

template<size_t D>
struct PiecewisePolynomial {
    std::vector<double> breakpoints;
    std::vector<Polynomial<D>> coefficients;
};

template<size_t D>
struct PolynomialPiece {
    Interval interval;
    Polynomial<D> polynomial;


};

template<size_t D>
void propagate_cost(const PiecewisePolynomial<D>& bottom);

/**
 * Given a bivariate polynomial function H(x,y) computes a univariate piecewise polynomial that returns for each Y
 * the minimum value bounded by the interval and left/right constraints (i.e. y maps to min_x H(x,y))
 *
 * @tparam D - Degree of the polynomials
 * @param h - Bivariate polynomial function H(x,y)
 * @param interval - Domain of Y
 * @param left_constraints - Constraints of shape a*y+b >= 0
 * @param right_constraints - Constraints of shape a*y+b <= 0
 * @return
 */
template<size_t D>
PiecewisePolynomial<D> find_minimum(
    const BivariatePolynomial<D>& h,
    const Interval& interval,
    std::vector<Polynomial<1>> left_constraints,
    std::vector<Polynomial<1>> right_constraints
) {
    // Minimum must lie on:
    //  - Derivative (w.r.t. x?) = 0 line (bounded by linear constraints)
    //  - Left/right boundary of the valid area (bounded by linear constraints)
    //
    // We take the cost function over these edges. This gives us a bunch of polynomial pieces over y.
    // The final result is the lower envelope of all theses pieces, which is a piecewise polynomial.

    // TODO: Add edges for critical lines (i.e. where partial derivative equals 0) of h.

    using Edge = PolynomialPiece<1>;
    std::vector<Edge> edges;

    std::vector<Polynomial<1>> lines;
    lines.insert(lines.end(), left_constraints.begin(), left_constraints.end());
    lines.insert(lines.end(), right_constraints.begin(), right_constraints.end());

    std::set<double> events;
    events.insert(interval.min);
    events.insert(interval.max);

    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 1; j < lines.size(); ++j) {
            auto intersections = find_intersections(lines[i], lines[j]);
            for (auto intersection : intersections) {
                if (interval.contains(intersection)) {
                    events.insert(intersection);
                }
            }
        }
    }

    bool prev_open = false;
    double left_start;
    Polynomial<1> left_current;
    double right_start;
    Polynomial<1> right_current;

    for (auto event : events) {
        auto compare = Polynomial<1>::CompareAt(event);
        std::sort(left_constraints.begin(), left_constraints.end(), compare);
        std::sort(right_constraints.begin(), right_constraints.end(), compare);

        auto left_max = left_constraints.back();
        auto right_min = right_constraints.front();

        bool is_open = !compare(right_min, left_max);

        if (prev_open) {
            // Close previously opened edges
            if (!is_open || left_max != left_current) {
                edges.push_back({{left_start, event}, left_current});
            }
            if (!is_open || right_min != right_current) {
                edges.push_back({{right_start, event}, right_current});
            }
        }

        if (is_open) {
            // Open new edges
            if (!prev_open || left_max != left_current) {
                left_start = event;
                left_current = left_max;
            }
            if (!prev_open || right_min != right_current) {
                right_start = event;
                right_current = right_min;
            }
        }

        if (prev_open && !is_open) {
            break;
        }

        prev_open = is_open;
    }

    for (auto edge : edges) {
        std::cout << edge.interval << ' ' << edge.polynomial << '\n';
    }


    std::vector<PolynomialPiece<D>> pieces;
    // TODO: Given a set of edges, compute the polynomial pieces

    return PiecewisePolynomial<D>();
}

template<size_t D>
PiecewisePolynomial<D> lower_envelope(std::vector<PolynomialPiece<D>>);
