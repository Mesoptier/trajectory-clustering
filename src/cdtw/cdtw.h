#pragma once

#include <vector>
#include <array>
#include <ostream>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include "Interval.h"
#include "Polynomial.h"
#include "BivariatePolynomial.h"

template<size_t D>
struct PiecewisePolynomial {
    std::vector<double> breakpoints;
    std::vector<Polynomial<D>> coefficients;
};

template<size_t D>
struct PolynomialPiece {
    Interval interval;
    Polynomial<D> polynomial;

    PolynomialPiece(const Interval& interval, const Polynomial<D>& polynomial) :
        interval(interval), polynomial(polynomial) {}

    friend std::ostream& operator<<(std::ostream& os, const PolynomialPiece& piece) {
        os << "{ " << piece.polynomial << ", " << piece.interval << " }";
        return os;
    }
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

    using Edge = PolynomialPiece<1>;
    std::vector<Edge> edges;

    std::vector<Polynomial<1>> lines;
    lines.insert(lines.end(), left_constraints.begin(), left_constraints.end());
    lines.insert(lines.end(), right_constraints.begin(), right_constraints.end());

    // TODO: Fix this for D != 2
    std::vector<Polynomial<1>> critical_lines = find_roots_y(h.partial_derivative_x());
    lines.insert(lines.end(), critical_lines.begin(), critical_lines.end());

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

    std::vector<bool> critical_open(critical_lines.size());
    std::vector<double> critical_start(critical_lines.size());

    for (auto event : events) {
        auto compare = Polynomial<1>::CompareAt(event);
        std::sort(left_constraints.begin(), left_constraints.end(), compare);
        std::sort(right_constraints.begin(), right_constraints.end(), compare);

        auto left_max = left_constraints.back();
        auto right_min = right_constraints.front();

        bool is_last_event = event == *events.rbegin();
        bool is_open = !compare(right_min, left_max) && !is_last_event;

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

        for (size_t i = 0; i < critical_lines.size(); ++i) {
            if (!critical_open[i]) {
                // Open new edges
                if (!compare(right_min, critical_lines[i]) && !compare(critical_lines[i], left_max)) {
                    critical_open[i] = true;
                    critical_start[i] = event;
                }
            } else {
                // Close previously opened edges
                if (is_last_event || compare(right_min, critical_lines[i]) || compare(critical_lines[i], left_max)) {
                    critical_open[i] = false;
                    edges.push_back({{critical_start[i], event}, critical_lines[i]});
                }
            }
        }

        if (prev_open && !is_open) {
            break;
        }

        prev_open = is_open;
    }

    std::vector<PolynomialPiece<D>> pieces;
    for (auto edge : edges) {
        const auto& hc = h.coefficients;
        const auto& ec = edge.polynomial.coefficients;

        // TODO: Fix this for D != 2
        pieces.emplace_back(edge.interval, Polynomial<2>({
            hc[0][0] + hc[1][0] * ec[0] + hc[2][0] * ec[0] * ec[0],
            hc[0][1] + hc[1][0] * ec[1] + hc[1][1] * ec[0] + 2 * hc[2][0] * ec[1] * ec[0],
            hc[0][2] + hc[1][1] * ec[1] + hc[2][0] * ec[1] * ec[1],
        }));
    }

    return lower_envelope(pieces);
}

template<size_t D>
PiecewisePolynomial<D> lower_envelope(const std::vector<PolynomialPiece<D>>& pieces) {
    for (auto piece : pieces) {
        std::cout << "Piecewise[{" << piece << "}, None],\n";
    }
}
