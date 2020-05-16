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
#include "PiecewisePolynomial.h"
#include "BivariatePolynomial.h"
#include "ConstrainedBivariatePolynomial.h"
#include "naive_lower_envelope.h"
#include "../IntegralFrechet/Cell.h"

template<size_t D>
void verify_minimum(const PiecewisePolynomial<D>& f, const std::vector<ConstrainedBivariatePolynomial<D>>& hs) {
    std::stringstream ss_h;
    for (const auto& h : hs) {
        ss_h << "Piecewise[{" << h << "}, None], ";
    }
    auto debug_h = ss_h.str();

    std::stringstream ss_f;
    ss_f << f;
    auto debug_f = ss_f.str();

    std::stringstream ss_m;
    ss_m << std::setprecision(30);
    ss_m << "{";

    bool brk = false;

    Interval f_interval = f.interval();
    size_t n = 100;
    for (size_t i = 1; i < n; ++i) {
        double y = f_interval.interpolate(i / static_cast<double>(n));
        double supposed_minimum = f(y);

        double real_minimum = std::numeric_limits<double>::infinity();
        for (const ConstrainedBivariatePolynomial<D>& h : hs) {
            if (h.y_interval.contains(y)) {
                auto slice = h.slice_at_y(y);
                real_minimum = std::min(real_minimum, slice.min_value());
            }
        }

        ss_m << "{" << y << ", " << real_minimum << "},";

        if (!std::isinf(real_minimum) && !approx_equal(supposed_minimum, real_minimum)) {
            brk = true;
        }
//        assert(approx_equal(supposed_minimum, real_minimum));
    }

    ss_m << "}";
    auto debug_m = ss_m.str();

    if (brk) {
        std::cout << "Plot3D[{" << debug_h << "} // Evaluate, {x,0,10}, {y,0,10}]\n";
        std::cout << "Show[Plot[{" << debug_f << "}, {x,0,10}], ListPlot[" << debug_m << "]]\n" << std::endl;
    }
//    assert(!brk);
}

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

    // TODO: Make this properly (instead of assuming 0y+0 is in left_constraints)
    // Early return for empty domain on x-axis
    for (const auto& rc : right_constraints) {
        if (rc.coefficients[0] == 0 && rc.coefficients[1] == 0) {
            return {};
        }
    }

    using Edge = PolynomialPiece<1>;
    std::vector<Edge> edges;

    std::vector<Polynomial<1>> lines;
    lines.insert(lines.end(), left_constraints.begin(), left_constraints.end());
    lines.insert(lines.end(), right_constraints.begin(), right_constraints.end());

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
        pieces.emplace_back(edge.interval, h.embed_x(edge.polynomial));
    }

    const PiecewisePolynomial<D> result = naive_lower_envelope(pieces);

    #ifndef NDEBUG
    if (!result.empty()) {
        // Checking whether the supposed minimum is actually the minimum
        const Interval result_interval = result.interval();
        const ConstrainedBivariatePolynomial<D> hc{
            h,
            interval,
            left_constraints,
            right_constraints,
        };

        size_t n = 100;
        for (size_t i = 0; i <= n; ++i) {
            double y = result_interval.interpolate(i / static_cast<double>(n));
            double supposed_minimum = result(y);
            auto slice = hc.slice_at_y(y);

            double real_minimum = slice.min_value();
            assert(approx_equal(supposed_minimum, real_minimum));
        }

    }
    #endif

    return result;
}

template<size_t D>
void fast_lower_envelope(PiecewisePolynomial<D>& left, const PiecewisePolynomial<D>& right) {
    // FIXME: Temporarily using naive_lower_envelope in order to reduce the chance of errors
    std::vector<PolynomialPiece<D>> pieces = left.pieces;
    pieces.insert(pieces.end(), right.pieces.begin(), right.pieces.end());
    left = naive_lower_envelope(pieces);
    return;

    // Early return for empty cases
    if (left.empty()) {
        left.pieces.assign(right.pieces.begin(), right.pieces.end());
        return;
    }
    if (right.empty()) {
        return;
    }

    // Verify that `left` is actually to the left of `right`
    assert(left.interval().min <= right.interval().min);
    assert(left.interval().max <= right.interval().max);
    // Verify that `left` and `right` have some overlap
    assert(left.interval().max >= right.interval().min);

    // Early return for trivial cases
    if (left.interval().max == right.interval().min) {
        assert(left.pieces.back().polynomial(left.interval().max)
                   == right.pieces.front().polynomial(right.interval().min));
        left.pieces.insert(left.pieces.end(), right.pieces.begin(), right.pieces.end());
        return;
    }
    if (right.interval().max == left.interval().max && left.pieces.back().polynomial(left.interval().max) < right.pieces.back().polynomial(right.interval().max)) {
        return;
    }

    // Find the rightmost piece in `right` which overlaps with the rightmost piece in `left`
    using PieceIterator = typename std::vector<PolynomialPiece<D>>::const_reverse_iterator;
    PieceIterator right_it = std::lower_bound(
        right.pieces.rbegin(),
        right.pieces.rend(),
        left.interval().max,
        [](const PolynomialPiece<D>& piece, double x) {
            return piece.interval.min >= x;
        }
    );


    // Find the point (X) where RIGHT is no longer lower than LEFT. This point must lie in the range where the intervals
    // of LEFT and RIGHT overlap. If no such point is found (i.e. either LEFT or RIGHT is always lower) we return the
    // lowest. Otherwise we split LEFT and RIGHT at X and fuse them together.

    double x = left.interval().max;

    while (!left.empty() && right_it != right.pieces.rend()) {
        const auto& left_back = left.pieces.back();

        assert(left_back.polynomial(x) >= right_it->polynomial(x));

        if (left_back.polynomial(x) == right_it->polynomial(x)) {
            return;
        }

        // Check intersections
        const auto intersections = find_intersections(left_back, *right_it);

        // Verify that there is at most one intersection between the two functions
        assert(intersections.size() <= 1);

        if (!intersections.empty()) {
            // Found an intersection between two pieces
            double root = intersections.front();

            // Shrink last piece of `left` or remove it if it would be empty
            if (left_back.interval.min == root) {
                left.pieces.pop_back();
            } else {
                left.pieces.back() = PolynomialPiece<D>({ left_back.interval.min, root }, left_back.polynomial);
            }

            // Add first piece from `right`
            left.pieces.push_back(PolynomialPiece<D>({ root, right_it->interval.max }, right_it->polynomial));
            left.pieces.insert(left.pieces.end(), right_it.base(), right.pieces.end());

            return;
        }

        if (left_back.interval.min == right_it->interval.min) {
            x = left_back.interval.min;
            left.pieces.erase(left.pieces.end() - 1);
            right_it++;
        } else if (left_back.interval.min < right_it->interval.min) {
            x = right_it->interval.min;
            right_it++;
        } else { // left_back.interval.min > right_it->interval.min
            x = left_back.interval.min;
            left.pieces.erase(left.pieces.end() - 1);
        }
    }

    // If left is empty, left and right were equal
    if (left.empty()) {
        left.pieces.assign(right.pieces.begin(), right.pieces.end());
        return;
    }

    throw std::logic_error("unresolved case");
}

template<size_t dimension, Norm image_norm, Norm param_norm>
class CDTW
{
public:
    static const size_t D =
        (dimension == 1 && image_norm == Norm::L1 && param_norm == Norm::L1) ? 2 :
        (dimension == 1 && image_norm == Norm::L2Squared && param_norm == Norm::L1) ? 3 :
        0;

    struct Entry
    {
        // Functions along x- and y-axis respectively
        PiecewisePolynomial<D> bottom;
        PiecewisePolynomial<D> left;
    };

private:
    const Curve& curve1;
    const Curve& curve2;

    size_t n;
    size_t m;

    // Dynamic program
    std::vector<std::vector<Entry>> in_functions;

    // Result
    PiecewisePolynomial<D> out_top;
    PiecewisePolynomial<D> out_right;

    // Internal methods
    template<class Iterator>
    PiecewisePolynomial<D> propagate(
        const std::vector<ConstrainedBivariatePolynomial<D>>& cell_costs,
        Iterator pieces_it,
        Iterator pieces_end
    ) const
    {
        std::vector<PolynomialPiece<D>> min_pieces;

        #ifndef NDEBUG
        std::vector<ConstrainedBivariatePolynomial<D>> total_costs;
        #endif

        PiecewisePolynomial<D> out_cost;
        while (pieces_it != pieces_end) {
            // piece_in_cost: cost of optimal path from origin to point on in-boundary
            const PolynomialPiece<D>& piece_in_cost = *pieces_it;

            // piece_out_cost: cost of optimal path from origin to point on out-boundary
            PiecewisePolynomial<D> piece_out_cost;

            // cell_cost: cost of optimal path from point on in-boundary to point on out-boundary
            for (const auto& cell_cost : cell_costs) {
                // total_cost: cost of optimal path from origin through point on in-boundary to point on out-boundary
                const auto total_cost = cell_cost.add_x(piece_in_cost);

                #ifndef NDEBUG
                total_costs.push_back(total_cost);
                #endif

                const PiecewisePolynomial<D> min_total_cost = find_minimum(
                    total_cost.f,
                    total_cost.y_interval,
                    total_cost.left_constraints,
                    total_cost.right_constraints
                );

                min_pieces.insert(min_pieces.end(), min_total_cost.pieces.begin(), min_total_cost.pieces.end());
            }

            ++pieces_it;
        }

        const auto result = naive_lower_envelope(min_pieces);

        #ifndef NDEBUG
        verify_minimum(result, total_costs);
        #endif

        return result;
    }
    PiecewisePolynomial<D> bottom_to_right(const PiecewisePolynomial<D>& in, const Cell& cell) const {
        // Walk through the in-pieces in reverse order. Because the first value in the out-function "originates from" the
        // last value in the in-function.
        return propagate(bottom_to_right_costs(cell), in.pieces.rbegin(), in.pieces.crend());
    }
    PiecewisePolynomial<D> bottom_to_top(const PiecewisePolynomial<D>& in, const Cell& cell) const {
        return propagate(bottom_to_top_costs(cell), in.pieces.begin(), in.pieces.cend());
    }

    // Internal methods specific to each case
    PiecewisePolynomial<D> base_bottom(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<D>> bottom_to_right_costs(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<D>> bottom_to_top_costs(const Cell& cell) const;

public:
    CDTW(const Curve& curve1, const Curve& curve2);

    double cost() const {
        const PiecewisePolynomial<D> func = in_functions[in_functions.size() - 2].back().bottom;
        return func.pieces.back().polynomial(func.interval().max);
        return 0;
    }

    void print_complexity() const {
        size_t max_complexity = 0;
        size_t total_complexity = 0;
        size_t num_functions = 0;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                size_t complexity_left = in_functions[i][j].left.pieces.size();
                size_t complexity_bottom = in_functions[i][j].bottom.pieces.size();

                max_complexity = std::max({ max_complexity, complexity_left, complexity_bottom });
                total_complexity += complexity_left + complexity_bottom;
                num_functions += (complexity_left != 0) + (complexity_bottom != 0);
            }
        }

        std::cout << "Max complexity: " << max_complexity << "\n";
        std::cout << "Total complexity: " << total_complexity << "\n";
        std::cout << "Average complexity: " << ((double) total_complexity / num_functions) << "\n";
    }

    void output_visualization_data() const {
        //
        // Visualization
        //

        std::ofstream file("data/out/exact-cdtw.txt");

        file << "Show[";
        for (size_t i = 0; i < n; ++i) {
            const double x = curve1.curve_length(i);

            for (size_t j = 0; j < m; ++j) {
                const double y = curve2.curve_length(j);

                const auto& cell_bottom = in_functions[i][j].bottom;
                if (!cell_bottom.empty()) {
                    file << "ParametricPlot3D[";
                    file << "{(x+" << x << "), (" << y << "), " << cell_bottom << "}";
                    file << ", {x," << cell_bottom.interval().min << "," << cell_bottom.interval().max << "}";
                    file << "],";
                }

                const auto& cell_left = in_functions[i][j].left;
                if (!cell_left.empty()) {
                    file << "ParametricPlot3D[";
                    file << "{(" << x << "), (x+" << y << "), " << cell_left << "}";
                    file << ",{x," << cell_left.interval().min << "," << cell_left.interval().max << "}";
                    file << "],";
                }
            }
        }
        file << "BoxRatios -> {Automatic, Automatic, 10}]" << std::endl;

        file.close();
    }

    const std::vector<std::vector<Entry>>& get_functions() const {
        return in_functions;
    }
};

// TODO: Remove temporary checks to find problematic cases for Sampson's proof.
template<size_t D>
void check_function(const PiecewisePolynomial<D>& f) {
    for (size_t i = 1; i < f.pieces.size(); ++i) {
        const PolynomialPiece<D> f1 = f.pieces.at(i - 1);
        const PolynomialPiece<D> f2 = f.pieces.at(i);
        const Polynomial<D> f1d = f1.polynomial.derivative();
        const Polynomial<D> f2d = f2.polynomial.derivative();

        if (f1d(f1.interval.max) < 0) {
            if (f2d(f2.interval.min) > 0) {
                std::cout << "negative to positive\n";
            }

            if (!approx_equal(f1d(f1.interval.max), f2d(f2.interval.min))) {
                std::cout << "not continuous after negative: Piecewise[{" << f1 << ", " << f2 << "}, None]\n";
            }
        }
    }
}

template<size_t dimension, Norm image_norm, Norm param_norm>
CDTW<dimension, image_norm, param_norm>::CDTW(const Curve& curve1, const Curve& curve2) :
    curve1(curve1), curve2(curve2),
    n(curve1.size()), m(curve2.size()),
    in_functions(n, std::vector<Entry>(m))
{
    // TODO: Simplify cell structs?

    // Initialize bottom base in_functions
    {
        double c = 0;
        for (size_t i = 0; i + 1 < n; ++i) {
            const Cell cell(curve1[i], curve2[0], curve1[i + 1], curve2[1]);
            const PiecewisePolynomial<D> f = base_bottom(cell) + c;
            in_functions[i][0].bottom = f;
            c = f.pieces.back().polynomial(f.pieces.back().interval.max); // TODO: Create and use back_value() instead
        }
    }

    // Initialize left base in_functions
    {
        double c = 0;
        for (size_t i = 0; i + 1 < m; ++i) {
            // Note: this cell is transposed
            const Cell cell(curve2[i], curve1[0], curve2[i + 1], curve1[1]);
            const PiecewisePolynomial<D> f = base_bottom(cell) + c;
            in_functions[0][i].left = f;
            c = f.pieces.back().polynomial(f.pieces.back().interval.max); // TODO: Create and use back_value() instead
        }
    }

    // Dynamic program
    for (size_t i = 0; i + 1 < n; ++i) {
        for (size_t j = 0; j + 1 < m; ++j) {
            // In-functions along bottom/left are given
            // Compute out-functions along top/right
            // Depending on current row/column:
            //  - set out-functions as in-function for next cell
            //  - or append to out_top/out_right

            // In-functions
            PiecewisePolynomial<D> bottom_in = in_functions[i][j].bottom;
            PiecewisePolynomial<D> left_in = in_functions[i][j].left;

            assert(!bottom_in.empty());
            assert(!left_in.empty());

            // Cell and transposed cell
            const Cell cell(curve1[i], curve2[j], curve1[i + 1], curve2[j + 1]);
            const Cell cell_t(curve2[j], curve1[i], curve2[j + 1], curve1[i + 1]);

            // Compute right
            PiecewisePolynomial<D> right_out = bottom_to_right(bottom_in, cell);
            fast_lower_envelope(right_out, bottom_to_top(left_in, cell_t));
            in_functions[i + 1][j].left = right_out;

            // Compute top
            PiecewisePolynomial<D> top_out = bottom_to_right(left_in, cell_t);
            fast_lower_envelope(top_out, bottom_to_top(bottom_in, cell));
            in_functions[i][j + 1].bottom = top_out;

//            // TODO: Remove temporary checks to find problematic cases for Sampson's proof.
//            check_function(right_out);
//            check_function(top_out);
        }
    }
}
