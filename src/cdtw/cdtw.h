#pragma once

#include <vector>
#include <array>
#include <ostream>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include <iomanip>

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

    //     if (!std::isinf(real_minimum) && !approx_equal(supposed_minimum, real_minimum)) {
    //         brk = true;
    //     }
    //    assert(approx_equal(supposed_minimum, real_minimum));
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
 * @param interval_y - Domain of Y
 * @param left_constraints - Constraints of shape a*y+b >= 0
 * @param right_constraints - Constraints of shape a*y+b <= 0
 * @return
 */
template<size_t D>
PiecewisePolynomial<D> find_minimum(
    const BivariatePolynomial<D>& h,
    const Interval& interval_y,
    std::vector<Polynomial<1>> left_constraints,
    std::vector<Polynomial<1>> right_constraints,
    const ConstrainedBivariatePolynomial<D>& f
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

    for (auto line: critical_lines) {
        double slope = line.coefficients[1];
        // std::cout << slope << std::endl;
    }

    std::set<double> events;
    events.insert(interval_y.min);
    events.insert(interval_y.max);

    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 1; j < lines.size(); ++j) {
            auto intersections = find_intersections(lines[i], lines[j]);
            for (auto intersection : intersections) {
                if (interval_y.contains(intersection)) {
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
        auto new_piece = PolynomialPiece<D>(edge.interval, h.embed_x(edge.polynomial));

        new_piece.history = History(
            f.path_type, f.boundaries, f.axis, edge.polynomial
        );

        auto pt = new_piece.history.path_type;

        if (pt != AU && pt != UA && pt != UAU && pt != AAU && pt != AAA && pt != UAA) {
            std::cout << "this is odd\n";
        }

        pieces.emplace_back(new_piece);
        // pieces.emplace_back(edge.interval, h.embed_x(edge.polynomial));
    }


    // std::cout << "number of pieces: " << pieces.size() << std::endl;

    return naive_lower_envelope(pieces);
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
        (image_norm == Norm::L1 && param_norm == Norm::L1) ? 2 :
        (image_norm == Norm::L2Squared && param_norm == Norm::L1) ? 3 :
        0;

    struct Entry
    {
        // Functions along x- and y-axis respectively
        PiecewisePolynomial<D> bottom;
        PiecewisePolynomial<D> left;
    };

private:
    Curve& curve1;
    Curve& curve2;

    size_t n;
    size_t m;

    double cdtw_cost;

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

        auto pieces = std::vector<PolynomialPiece<D>>();
        auto new_it = pieces_it;
        while (new_it != pieces_end) {
            pieces.push_back(*new_it);
            ++new_it;
        }
        int counter = 0;
        PiecewisePolynomial<D> out_cost;
        while (pieces_it != pieces_end) {
            // piece_in_cost: cost of optimal path from origin to point on in-boundary
            const PolynomialPiece<D>& piece_in_cost = *pieces_it;
            counter++;
            // piece_out_cost: cost of optimal path from origin to point on out-boundary
            PiecewisePolynomial<D> piece_out_cost;

            // if (approx_equal(cell_costs[0].f.coefficients[0][0], 0.0194395422591176) && approx_equal(cell_costs[0].f.coefficients[0][1], 0.18967521022380202)) {
            //     std::cout << "hi\n";        
            // }

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
                    total_cost.right_constraints,
                    cell_cost
                );

                auto single_variable = total_cost.slice_at_y(1.4142630607139537);

                auto single_variable_cell_cost = cell_cost.slice_at_y(1.4142630607139537);

                // if (min_total_cost.pieces.size() > 0) {
                //     auto poly = min_total_cost.pieces.back();
                //     if (approx_equal(poly.polynomial(poly.interval.max) / 100000, 8.21263 / 100000))
                //         std::cout << "hello\n";
                //     else
                //         std::cout << "value: " << poly.polynomial(poly.interval.max) << " const: " << poly.polynomial.coefficients[0] << std::endl;
                // }

                

                min_pieces.insert(min_pieces.end(), min_total_cost.pieces.begin(), min_total_cost.pieces.end());
            }

            ++pieces_it;
        }
        // std::cout << "count: " << counter << std::endl;
        const auto result = naive_lower_envelope(min_pieces);


        // if (approx_equal(cell_cost.f.coefficients[0][0], 1.4141994337255621) 
        //         && approx_equal(cell_cost.f.coefficients[2][0], 7.0710501343442271e-06))
        //             std::cout << "hi\n";
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

    void get_path_segs(PathType path_type, Boundaries boundaries, double in, double out,
    std::vector<std::vector<double>>& output, Curve& c1, Curve& c2, int i, int j, Line axis, bool left) {

        std::vector<double> e1;
        std::vector<double> e2;
        std::vector<double> e3;
        std::vector<double> e4;

        double width; 
        if (left && j > 0)
            width = c2.curve_length(j-1, j);
        else if (j < c2.size()-1)
            width = c2.curve_length(j, j+1);
        else
            width = 10000000;

        double height;
        if (left)
            height = c1.curve_length(i, i+1);
        else if (i > 0)
            height = c1.curve_length(i-1, i);
        else
            height = 10000000;

        std::vector<double> left_to_axis_a = {0, in, axis.getX(in), in};
        std::vector<double> bottom_to_axis_u = {in, 0, in, axis.getY(in)};
        std::vector<double> axis_to_right_a = {axis.getX(out), out, width, out};
        std::vector<double> axis_to_top_u = {out, axis.getY(out), out, height};

        std::vector<double> left_to_axis_u = {0, in, 0, axis.getY(0)};
        std::vector<double> bottom_to_axis_a = {in, 0, axis.getX(0), 0};
        std::vector<double> axis_to_right_u = {width, axis.getY(width), width, out};
        std::vector<double> axis_to_top_a = {axis.getX(height), height, out, height};

        std::vector<double> along_bt;// = {in, axis.getY(in), out, axis.getY(out)};
        std::vector<double> along_br;// = {in, axis.getY(in), axis.getX(out), out};
        std::vector<double> along_lr;// = {axis.getX(in), in, axis.getX(out), out};
        std::vector<double> along_lt;// = {axis.getX(in), in, out, axis.getY(out)};

        switch (path_type) {
            case (UAU): {
                along_bt = {in, axis.getY(in), out, axis.getY(out)};
                along_br = {in, axis.getY(in), axis.getX(out), out};
                along_lr = {0, axis.getY(0), width, axis.getY(width)};
                along_lt = {0, axis.getY(0), out, axis.getY(out)};
                break;
            }
            case (UAA): {
                along_bt = {in, axis.getY(in), axis.getX(height), height};
                along_br = {in, axis.getY(in), axis.getX(out), out};
                along_lr = {0, axis.getY(0), axis.getX(out), out};
                along_lt = {0, axis.getY(0), axis.getX(height), height};
                break;
            }
            case (AAA): {
                along_bt = {axis.getX(0), 0, axis.getX(height), height};
                along_br = {axis.getX(0), 0, width, axis.getY(width)};
                along_lr = {axis.getX(in), in, axis.getX(out), out};
                along_lt = {axis.getX(in), in, axis.getX(height), height};
                break;
            }
            case (AAU): {
                along_bt = {axis.getX(0), 0, out, axis.getY(out)};
                along_br = {axis.getX(0), 0, axis.getX(out), out};
                along_lr = {axis.getX(in), in, width, axis.getY(width)};
                along_lt = {axis.getX(in), in, out, axis.getY(out)};
                break;
            }
        }


        if (path_type == UA) {

            if (boundaries == BR) {

                e1 = {width, out, in, out};
                e2 = {in, out, in, 0};

            } else if (boundaries == BT) {

                e1 = {out, height, in, height};
                e2 = {in, height, in, 0};

            } else if (boundaries == LT) {

                e1 = {out, height, 0, height};
                e2 = {0, height, 0, in};

            } else if (boundaries == LR) {

                e1 = {width, out, 0, out};
                e2 = {0, in, 0, in};

            }

        } else if (path_type == AU) {

            if (boundaries == BR) {

                e1 = {width, out, width, 0};
                e2 = {width, 0, in, 0};

            } else if (boundaries == BT) {

                e1 = {out, height, out, 0};
                e2 = {out, 0, in, 0};

            } else if (boundaries == LT) {

                e1 = {out, height, 0, height};
                e2 = {0, height, 0, in};

            } else if (boundaries == LR) {

                e1 = {width, out, 0, out};
                e2 = {0, out, 0, in};

            }

        } else {
            
            if (boundaries == BR) {
                
                if (path_type == UAU) {

                    e1 = bottom_to_axis_u;
                    e2 = along_br;
                    e3 = axis_to_right_u;

                } else if (path_type == UAA) {

                    e1 = bottom_to_axis_u;
                    e2 = along_br;
                    e3 = axis_to_right_a;

                } else if (path_type == AAA) {

                    e1 = bottom_to_axis_a;
                    e2 = along_br;
                    e3 = axis_to_right_a;

                } else if (path_type == AAU) {

                    e1 = bottom_to_axis_a;
                    e2 = along_br;
                    e3 = axis_to_right_u;

                }

            } else if (boundaries == BT) {

                if (path_type == UAU) {

                    e1 = bottom_to_axis_u;
                    e2 = along_bt;
                    e3 = axis_to_top_u;

                } else if (path_type == UAA) {

                    e1 = bottom_to_axis_u;
                    e2 = along_bt;
                    e3 = axis_to_top_a;

                } else if (path_type == AAA) {
                    
                    e1 = bottom_to_axis_a;
                    e2 = along_bt;
                    e3 = axis_to_top_a;

                } else if (path_type == AAU) {

                    e1 = bottom_to_axis_a;
                    e2 = along_bt;
                    e3 = axis_to_top_u;

                }

            } else if (boundaries == LR) {

                if (path_type == UAU) {

                    e1 = left_to_axis_u;
                    e2 = along_lr;
                    e3 = axis_to_right_u;

                } else if (path_type == UAA) {

                    e1 = left_to_axis_u;
                    e2 = along_lr;
                    e3 = axis_to_right_a;

                } else if (path_type == AAA) {
                    
                    e1 = left_to_axis_a;
                    e2 = along_lr;
                    e3 = axis_to_right_a;

                } else if (path_type == AAU) {

                    e1 = left_to_axis_a;
                    e2 = along_lr;
                    e3 = axis_to_right_u;

                }

            } else if (boundaries == LT) {

                if (path_type == UAU) {

                    e1 = left_to_axis_u;
                    e2 = along_lt;
                    e3 = axis_to_top_u;

                } else if (path_type == UAA) {

                    e1 = left_to_axis_u;
                    e2 = along_lt;
                    e3 = axis_to_top_a;

                } else if (path_type == AAA) {
                    
                    e1 = left_to_axis_a;
                    e2 = along_lt;
                    e3 = axis_to_top_a;

                } else if (path_type == AAU) {

                    e1 = left_to_axis_a;
                    e2 = along_lt;
                    e3 = axis_to_top_u;

                }

            }

        }

        if (((boundaries == BR) || (boundaries == BT)) && i == 0) {
            e4 = {0, 0, in, 0};
        } else if (((boundaries == LR) || (boundaries == LT)) && j == 0) {
            e4 = {0, 0, 0, in};
        }

        int count = 0;

        if (e1.size() > 0) {
            output.push_back(e1);
            ++count;
        }
        if (e2.size() > 0){
            output.push_back(e2);
            ++count;
        }
        if (e3.size() > 0) {
            output.push_back(e3);
            ++count;
        }
        if (e4.size() > 0) {
            output.push_back(e4);
            ++count;
        }

        double x_translation;
        double y_translation;

        for (int k = 0; k < count; ++k) {
            if (i > 0) {
                y_translation = left ? c1.curve_length(0, i) : c1.curve_length(0, i-1);
                output[output.size()-1-k][1] = output[output.size()-1-k][1] + y_translation;
                output[output.size()-1-k][3] = output[output.size()-1-k][3] + y_translation;
            }

            if (j > 0) {
                x_translation = left ? c2.curve_length(0, j-1) : c2.curve_length(0, j);
                output[output.size()-1-k][0] = output[output.size()-1-k][0] + x_translation;
                output[output.size()-1-k][2] = output[output.size()-1-k][2] + x_translation;
            }
        }

        std::cout << "";

    }

    void write_heat_map(Curve& new_curve1, Curve& new_curve2) {
        std::ofstream output("heat_map.dat");

        output << std::fixed << std::setprecision(10);

        auto len1 = new_curve1.curve_length();
        auto len2 = new_curve2.curve_length();

        int grid_size = 100;

        double x_step = len1 / grid_size;
        double y_step = len2 / grid_size;

        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j) {
                CPoint cp1 = new_curve1.get_cpoint_after(x_step * i);
                CPoint cp2 = new_curve2.get_cpoint_after(y_step * j);

                PointID p1_id = cp1.getPoint();
                distance_t p1_frac = cp1.getFraction();
                Point p1 = new_curve1[p1_id] + (new_curve1[p1_id+1] - new_curve1[p1_id])*p1_frac;

                std::cout << p1_id << " " << p1_frac << "\n";
                std::cout << p1 << "\n";

                PointID p2_id = cp2.getPoint();
                distance_t p2_frac = cp2.getFraction();
                Point p2 = new_curve2[p2_id] + (new_curve2[p2_id+1] - new_curve2[p2_id])*p2_frac;


                // Point p1 = new_curve1[cp1.getPoint()];
                // Point p2 = new_curve2[cp2.getPoint()];
                double height = abs(p1.x - p2.x) + abs(p1.y - p2.y);


                output << height;
                if (j < grid_size - 1)
                    output << " ";
                else
                    output << "\n";
            }


        output.close();
    }

    void compute_warping_path(Curve& new_curve1, Curve& new_curve2, 
    std::vector<std::vector<double>>& output, int i, int j, double y, bool left) {

        if (i + j == 0)
            return;

        PolynomialPiece<2> piece = PolynomialPiece<2>({0, 1}, Polynomial<2>({{0, 0, 0}}));

        if (left) {
            piece = in_functions[i][j].left.get_piece(y);
        } else {
            piece = in_functions[i][j].bottom.get_piece(y);
        }

        auto boundaries = piece.history.boundaries;
        auto path_type = piece.history.path_type;
        auto axis = piece.history.axis;
        double x = piece.history.x_poly(y);

        std::cout << "i: " << i << " j: " << j << std::endl;
        get_path_segs(path_type, boundaries, x, y, output, new_curve1, new_curve2, i, j, axis, left);

        

        switch(boundaries) {
            case BR:
                assert(left);
                if (j == 0)
                    break;
                else
                    compute_warping_path(new_curve1, new_curve2, output, i, j-1, x, false);
                break;
            case LR:
                assert(left);
                if (j == 0)
                    break;
                else
                    compute_warping_path(new_curve1, new_curve2, output, i, j-1, x, true);
                break;
            case BT:
                assert(!left);
                if (i == 0)
                    break;
                else
                    compute_warping_path(new_curve1, new_curve2, output, i-1, j, x, false);
                break;
            case LT:
                assert(!left);
                if (i == 0)
                    break;
                else
                    compute_warping_path(new_curve1, new_curve2, output, i-1, j, x, true);
                break;
        }

        io::write_path("warping_path.txt", output, new_curve1.size(), new_curve2.size()); 
    }

    void process_cell(int i, int j, Curve& curve1, Curve& curve2) {
        
        auto cell_t = Cell(curve1[i], curve2[j], curve1[i+1], curve2[j+1]);
        const Cell cell(curve2[j], curve1[i], curve2[j+1], curve1[i+1]);
        
        if (i == 0 && j == 0) {
            auto bottom_poly = base_bottom(cell);
            auto left_poly = base_bottom(cell_t);

            for (auto piece: left_poly.pieces)
                piece.history.transpose();

            in_functions[0][0].bottom = bottom_poly;
            in_functions[0][0].left = left_poly;
        
        } 

        auto right_out_functions = std::vector<PolynomialPiece<D>>();
        auto top_out_functions = std::vector<PolynomialPiece<D>>();

        if (i > 0 || j == 0) {
            // std::cout << in_functions[i][j].bottom.pieces.size();

            auto b2r = bottom_to_right(in_functions[i][j].bottom, cell);
            for (auto piece: b2r.pieces)
                right_out_functions.push_back(piece);

            auto b2t = bottom_to_top(in_functions[i][j].bottom, cell);
            for (auto piece: b2t.pieces)
                top_out_functions.push_back(piece);

        }

        if (j > 0 || i == 0) {
            // std::cout << in_functions[i][j].left.pieces.size();
            if (i == 0 && j == 20)
                std::cout << "hello\n";
            auto l2r = bottom_to_top(in_functions[i][j].left, cell_t);
            for (auto piece: l2r.pieces) {
                right_out_functions.push_back(piece);
                right_out_functions.back().history.transpose();
            }


            auto l2t = bottom_to_right(in_functions[i][j].left, cell_t);
            for (auto piece: l2t.pieces) {
                top_out_functions.push_back(piece);
                top_out_functions.back().history.transpose();
            }
        }


        if (j < in_functions[i].size()-1)
            in_functions[i][j+1].left = naive_lower_envelope(right_out_functions);
        if (i < in_functions.size()-1)
            in_functions[i+1][j].bottom = naive_lower_envelope(top_out_functions);
    }

    // Internal methods specific to each case
    PiecewisePolynomial<D> base_bottom(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<D>> bottom_to_right_costs(const Cell& cell) const;
    std::vector<ConstrainedBivariatePolynomial<D>> bottom_to_top_costs(const Cell& cell) const;

public:
    CDTW(Curve& _curve1, Curve& _curve2);

    double cost() const {
        // return cdtw_cost;
        const PiecewisePolynomial<D> func = in_functions[in_functions.size() - 2].back().left;
        // if (func.pieces.back().polynomial(func.interval().max) < 0)
        //     std::cout << "hi\n";
        return func.pieces.back().polynomial(func.interval().max);
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
        file << "BoxRatios -> {Automatic, Automatic, 10}, PlotRange -> {Automatic, Automatic, {0, 1500}}]" << std::endl;

        file.close();
    }

    const std::vector<std::vector<Entry>>& get_functions() const {
        return in_functions;
    }

    void check_table() {
        for (int i = 1; i < in_functions.size(); ++i) {
            for (int j = 1; j < in_functions[0].size(); ++j) {
                std::cout << i << ", " << j << std::endl;
                auto in_left = in_functions[i][j].left.pieces;
                auto in_bottom = in_functions[i][j].bottom.pieces;

                std::cout << "bottom left corner: \n" <<
                "left: " << in_left[0].polynomial(0) << std::endl <<
                "bottom: " << in_bottom[0].polynomial(0) << std::endl;

                if (i > 0 || (i == 1 && j == 0)) {
                    std::cout << "lower cell, left top: " << in_functions[i-1][j].left.pieces.back().polynomial(
                        in_functions[i-1][j].left.pieces.back().interval.max
                    ) << std::endl;
                }

                if (j > 1 || (j == 1 && i == 0)) {
                    std::cout << "left cell, bottom right: " << in_functions[i][j-1].bottom.pieces.back().polynomial(
                        in_functions[i][j-1].bottom.pieces.back().interval.max
                    ) << std::endl;
                }

                // if(!approx_equal(in_left[0].polynomial(0)/100, in_bottom[0].polynomial(0)/100))
                //     assert(false);
            }
        }

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

namespace {
    enum class Side {
        LEFT,
        BOTTOM,
    };

    std::ostream& operator<<(std::ostream& os, const Side& side) {
        switch (side) {
            case Side::LEFT:
                os << "LEFT";
                break;
            case Side::BOTTOM:
                os << "BOTTOM";
                break;
        }
        return os;
    }

    struct Coord {
        size_t i1;
        size_t i2;
        Side side;

        Coord top() const {
            return {i1, i2 + 1, Side::BOTTOM};
        }
        Coord right() const {
            return {i1 + 1, i2, Side::LEFT};
        }

        struct hash
        {
            auto operator()(const Coord& x) const {
                size_t seed = 0;
                hash_combine(seed, x.i1);
                hash_combine(seed, x.i2);
                hash_combine(seed, x.side);
                return seed;
            }
        };

        bool operator==(const Coord& rhs) const {
            return i1 == rhs.i1 &&
                i2 == rhs.i2 &&
                side == rhs.side;
        }

        bool operator<(const Coord& rhs) const {
            if (i1 < rhs.i1)
                return true;
            if (rhs.i1 < i1)
                return false;
            if (i2 < rhs.i2)
                return true;
            if (rhs.i2 < i2)
                return false;
            return side < rhs.side;
        }

        friend std::ostream& operator<<(std::ostream& os, const Coord& coord) {
            os << "Coord{" << coord.i1 << ", " << coord.i2 << ", " << coord.side << "}";
            return os;
        }
    };

    struct PQNode {
        distance_t min_cost;
        Coord coord;

        bool operator<(const PQNode& rhs) const {
            if (min_cost < rhs.min_cost)
                return true;
            if (rhs.min_cost < min_cost)
                return false;
            return coord < rhs.coord;
        }
        bool operator>(const PQNode& rhs) const {
            return rhs < *this;
        }
        bool operator<=(const PQNode& rhs) const {
            return !(rhs < *this);
        }
        bool operator>=(const PQNode& rhs) const {
            return !(*this < rhs);
        }

        friend std::ostream& operator<<(std::ostream& os, const PQNode& node) {
            os << "PQNode{" << node.min_cost << ", " << node.coord << "}";
            return os;
        }
    };
}

template<size_t dimension, Norm image_norm, Norm param_norm>
CDTW<dimension, image_norm, param_norm>::CDTW(Curve& _curve1, Curve& _curve2) :
    curve1(_curve1), curve2(_curve2),
    n(curve1.size()), m(curve2.size()),
    in_functions(n, std::vector<Entry>(m))
{
    int new_n;
    int new_m;

    std::vector<CPoint> new_c_points_1 = std::vector<CPoint>();
    std::vector<CPoint> new_c_points_2 = std::vector<CPoint>();

    for (int i = 0; i < curve1.size(); ++i) {
        new_c_points_1.push_back(CPoint(i, 0));
    }

    for (int i = 0; i < curve2.size(); ++i) {
        new_c_points_2.push_back(CPoint(i, 0));
    }


    for (int i = 0; i < curve1.size()-1; ++i) {
        auto s1 = curve1[i];
        auto t1 = curve1[i+1];

        const auto l1 = Line::fromTwoPoints(s1, t1);

        for (int j = 0; j < curve2.size()-1; ++j) {

            auto s2 = curve2[j];
            auto t2 = curve2[j+1];

            const auto l2 = Line::fromTwoPoints(s2, t2);

            if (!isParallel(l1, l2)) {
                const auto p = intersect(l1, l2);
                Point mid = p;

                distance_t c1_dist = s1.dist(mid);
                distance_t c2_dist = s2.dist(mid);

                if ((mid.x - s1.x)*(mid.x - t1.x) < 0
                && (mid.x - s2.x)*(mid.x - t2.x) < 0) {
                    std::cout << c1_dist/(s1.dist(t1)) << std::endl;
                    std::cout << c2_dist/(s2.dist(t2)) << std::endl;
                    // std::cout << s1 << std::endl;
                    // std::cout << t1 << std::endl;
                    // std::cout << s2 << std::endl;
                    // std::cout << t2 << std::endl;
                    // std::cout << mid << std::endl;
                    new_c_points_1.push_back(CPoint(i, c1_dist/(s1.dist(t1))));
                    new_c_points_2.push_back(CPoint(j, c2_dist/(s2.dist(t2))));
                }
            }
        }
    }

    std::sort(new_c_points_1.begin(), new_c_points_1.end());
    std::sort(new_c_points_2.begin(), new_c_points_2.end());

    Points new_points_1 = Points();

    for (auto c: new_c_points_1)
        new_points_1.push_back(curve1.interpolate_at(c));

    Points new_points_2 = Points();

    for (auto c: new_c_points_2)
        new_points_2.push_back(curve2.interpolate_at(c));

    Curve new_curve1 = Curve("", new_points_1);
    Curve new_curve2 = Curve("", new_points_2);

    n = new_curve1.size();
    m = new_curve2.size();

    for (int i = 0; i < n-1; ++i)
        for (int j = 0; j < m-1; ++j) {
            if (i == 1 && j == 0)
                std::cout << "hi\n";
            process_cell(i, j, new_curve1, new_curve2);
        }

    auto output = std::vector<std::vector<double>>();
    write_heat_map(new_curve1, new_curve2);
    compute_warping_path(new_curve1, new_curve2, output, n-2, m-1, new_curve1.curve_length(n-2, n-1), true);
    io::export_points("c1.txt", new_curve1.get_points());
    io::export_points("c2.txt", new_curve2.get_points());
    std::cout << "segments: " << output.size() << std::endl;
    std::cout << "computed cdtw...\n";
}

// template<size_t dimension, Norm image_norm, Norm param_norm>
// CDTW<dimension, image_norm, param_norm>::CDTW(Curve& _curve1, Curve& _curve2) :
//     curve1(_curve1), curve2(_curve2),
//     n(curve1.size()), m(curve2.size()),
//     in_functions(n, std::vector<Entry>(m))
// {
//     struct Function {
//         PiecewisePolynomial<D> f;
//         distance_t min_cost;
//         bool open; // Whether this function is currently in the open_set
//     };

//     // Curve new_curve1("", curve1.get_points());
//     // Curve new_curve2("", curve2.get_points());

//     // if (dimension == 2 && image_norm == Norm::L1 && param_norm == Norm::L1) {

//         int new_n;
//         int new_m;

//         std::vector<CPoint> new_c_points_1 = std::vector<CPoint>();
//         std::vector<CPoint> new_c_points_2 = std::vector<CPoint>();

//         for (int i = 0; i < curve1.size(); ++i) {
//             new_c_points_1.push_back(CPoint(i, 0));
//         }

//         for (int i = 0; i < curve2.size(); ++i) {
//             new_c_points_2.push_back(CPoint(i, 0));
//         }


//         for (int i = 0; i < curve1.size()-1; ++i) {
//             auto s1 = curve1[i];
//             auto t1 = curve1[i+1];

//             const auto l1 = Line::fromTwoPoints(s1, t1);

//             for (int j = 0; j < curve2.size()-1; ++j) {

//                 auto s2 = curve2[j];
//                 auto t2 = curve2[j+1];

//                 const auto l2 = Line::fromTwoPoints(s2, t2);

//                 if (!isParallel(l1, l2)) {
//                     const auto p = intersect(l1, l2);
//                     Point mid = p;

//                     distance_t c1_dist = s1.dist(mid);
//                     distance_t c2_dist = s2.dist(mid);

//                     if ((mid.x - s1.x)*(mid.x - t1.x) < 0
//                     && (mid.x - s2.x)*(mid.x - t2.x) < 0) {
//                         std::cout << c1_dist/(s1.dist(t1)) << std::endl;
//                         std::cout << c2_dist/(s2.dist(t2)) << std::endl;
//                         // std::cout << s1 << std::endl;
//                         // std::cout << t1 << std::endl;
//                         // std::cout << s2 << std::endl;
//                         // std::cout << t2 << std::endl;
//                         // std::cout << mid << std::endl;
//                         new_c_points_1.push_back(CPoint(i, c1_dist/(s1.dist(t1))));
//                         new_c_points_2.push_back(CPoint(j, c2_dist/(s2.dist(t2))));
//                     }
//                 }
//             }
//         }

//         std::sort(new_c_points_1.begin(), new_c_points_1.end());
//         std::sort(new_c_points_2.begin(), new_c_points_2.end());
    
//         Points new_points_1 = Points();

//         for (auto c: new_c_points_1)
//             new_points_1.push_back(curve1.interpolate_at(c));

//         Points new_points_2 = Points();

//         for (auto c: new_c_points_2)
//             new_points_2.push_back(curve2.interpolate_at(c));

//         Curve new_curve1 = Curve("", new_points_1);
//         Curve new_curve2 = Curve("", new_points_2);
//     // }

//     std::unordered_map<Coord, Function, typename Coord::hash> functions;
//     std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> open_set;

//     {
//         const Cell cell(new_curve1[0], new_curve2[0], new_curve1[1], new_curve2[1]);

//         const Function f1{base_bottom(cell), 0, true};
//         const Coord c1{0, 0, Side::BOTTOM};
//         functions.insert_or_assign(c1, f1);
//         open_set.push({f1.f.min_value(), c1});
//     }

//     distance_t min_result_cost = std::numeric_limits<distance_t>::infinity();

//     size_t nodes_handled = 0;
//     size_t nodes_skipped = 0;
//     size_t nodes_opened = open_set.size();
//     size_t nodes_reopened = 0;

//     while (!open_set.empty()) {
//         const PQNode node = open_set.top();
//         open_set.pop();

//         std::cout << node.coord << std::endl;

// //        std::cout << min_result_cost << " " << node << '\n';
//         const Coord coord = node.coord;
//         Function& f = functions.at(coord);
//         f.open = false;

//         // We can no longer find a lower min_result_cost, so the search is finished
//         if (node.min_cost >= min_result_cost) {
//             break;
//         }

//         if (node.min_cost > f.min_cost) {
//             ++nodes_skipped;
//             continue;
//         }

//         if (coord.i1 == new_curve1.size() - 1 || coord.i2 == new_curve2.size() - 1) {
//             if (coord.i1 == new_curve1.size() - 2 || coord.i2 == new_curve2.size() - 2) {
//                 min_result_cost = std::min(min_result_cost, f.f.right_value());
//             }
//             ++nodes_skipped;
//             continue;
//         }

//         ++nodes_handled;

//         // Cell and transposed cell
//         const Cell cell(new_curve1[coord.i1], new_curve2[coord.i2], new_curve1[coord.i1 + 1], new_curve2[coord.i2 + 1]);
//         const Cell cell_t(new_curve2[coord.i2], new_curve1[coord.i1], new_curve2[coord.i2 + 1], new_curve1[coord.i1 + 1]);

//         for (size_t i = 0; i < 2; ++i) { // 0 -> right, 1 -> top
//             const PiecewisePolynomial<D> out_f = i == 0
//                 ? (coord.side == Side::BOTTOM ? bottom_to_right(f.f, cell) : bottom_to_top(f.f, cell_t))
//                 : (coord.side == Side::BOTTOM ? bottom_to_top(f.f, cell) : bottom_to_right(f.f, cell_t));
//             const Coord out_coord = i == 0
//                 ? coord.right()
//                 : coord.top();

            

//             distance_t min_cost = out_f.min_value();

//             auto it = functions.find(out_coord);
//             if (it == functions.end()) {
//                 functions.emplace(out_coord, Function{out_f, min_cost, true});
//                 open_set.push({min_cost, out_coord});
//                 ++nodes_opened;
//             } else {
//                 Function& prev_out = it->second;
//                 // TODO: Only re-open function if lower_envelope has changed
//                 fast_lower_envelope(prev_out.f, out_f);

//                 if (min_cost < prev_out.min_cost) {
//                     prev_out.open = true;
//                     open_set.push({min_cost, out_coord});
//                     ++nodes_reopened;
//                 } else { // min_cost == prev_out.min_cost
//                     if (!prev_out.open) {
//                         prev_out.open = true;
//                         open_set.push({min_cost, out_coord});
//                         ++nodes_reopened;
//                     }
//                 }
//             }
//         }
//     }


//     cdtw_cost = min_result_cost;
//     std::cout << min_result_cost << "\n";
//     std::cout << "nodes_handled: " << nodes_handled << "\n";
//     std::cout << "nodes_skipped: " << nodes_skipped << "\n";
//     std::cout << "nodes_opened: " << nodes_opened << "\n";
//     std::cout << "nodes_reopened: " << nodes_reopened << "\n";
//     std::cout << "nodes remaining: " << open_set.size() << "\n";
// }
