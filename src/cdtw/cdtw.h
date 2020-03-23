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
struct PolynomialPiece
{
    Interval interval;
    Polynomial<D> polynomial;

    PolynomialPiece(const Interval& interval, const Polynomial<D>& polynomial) :
        interval(interval), polynomial(polynomial) {}

    //
    // Arithmetic operators
    //

    PolynomialPiece<D>& operator+=(double c) {
        polynomial += c;
        return *this;
    }
    friend PolynomialPiece<D> operator+(PolynomialPiece<D> f, double c) {
        f += c;
        return f;
    }

    //
    // Stream output operator
    //

    friend std::ostream& operator<<(std::ostream& os, const PolynomialPiece& piece) {
        os << "{ " << piece.polynomial << ", " << piece.interval << " }";
        return os;
    }
};

template<size_t D>
std::vector<double> find_intersections(const PolynomialPiece<D>& f, const PolynomialPiece<D>& g) {
    std::vector<double> result;
    const std::vector<double> intersections = find_intersections(f.polynomial, g.polynomial);
    for (auto x : intersections) {
        if (f.interval.contains(x) && g.interval.contains(x)) {
            result.push_back(x);
        }
    }
    return result;
}

template<size_t D>
struct PiecewisePolynomial
{
    std::vector<PolynomialPiece<D>> pieces;

    PiecewisePolynomial() = default;
    PiecewisePolynomial(const std::vector<PolynomialPiece<D>>& pieces) : pieces(pieces) {
        for (size_t i = 1; i < pieces.size(); ++i) {
            const PolynomialPiece<D>& p1 = pieces[i - 1];
            const PolynomialPiece<D>& p2 = pieces[i];
            // Verify that piece intervals line up
            assert(p1.interval.max == p2.interval.min);
            // Verify that pieces connect
            assert(p1.polynomial(p1.interval.max) == p2.polynomial(p2.interval.min));
        }
    }

    bool empty() const {
        return pieces.empty();
    }

    Interval interval() const {
        return {
            pieces.front().interval.min,
            pieces.back().interval.max,
        };
    }

    //
    // Arithmetic operators
    //

    PiecewisePolynomial<D>& operator+=(double c) {
        for (PolynomialPiece<D>& piece : pieces) {
            piece += c;
        }
        return *this;
    }
    friend PiecewisePolynomial<D> operator+(PiecewisePolynomial<D> f, double c) {
        f += c;
        return f;
    }

    //
    // Stream output operator
    //

    friend std::ostream& operator<<(std::ostream& os, const PiecewisePolynomial& f) {
        os << "Piecewise[{";
        for (size_t i = 0; i < f.pieces.size(); ++i) {
            if (i != 0) os << ",";
            os << f.pieces[i];
        }
        os << "}, None]";
        return os;
    }
};

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

    return naive_lower_envelope(pieces);
}

/**
 * Given a set of polynomial pieces, computes the set of polynomial pieces that form the lower envelope without any
 * additional knowledge about the pieces (such as the crossing lemma).
 *
 * @tparam D
 * @param pieces
 * @return
 */
template<size_t D>
PiecewisePolynomial<D> naive_lower_envelope(const std::vector<PolynomialPiece<D>>& pieces) {
//    for (auto piece : pieces) {
//        std::cout << "Piecewise[{" << piece << "}, None],\n";
//    }

    // State:
    //  - Open pieces sorted by value (and derivatives) at x
    // Events:
    //  - Piece opens/closes: update open pieces list
    //  - Intersection between open pieces: re-order open piece list

    // PIECE_OPEN: id = opened
    // PIECE_CLOSE: id = closed
    // PIECE_INTERSECT: id = lowest of intersecting pieces after intersection

    using PieceID = size_t;
    constexpr PieceID INVALID_PIECE_ID = std::numeric_limits<size_t>::max();
    enum class EventType
    {
        OPEN,
        CLOSE,
        SWAP,
    };
    struct Event
    {
        double x;
        EventType type;
        PieceID id;
    };
    struct CompareEvent
    {
        const std::vector<PolynomialPiece<D>>& pieces;
        explicit CompareEvent(const std::vector<PolynomialPiece<D>>& pieces) : pieces(pieces) {}
        bool operator()(const Event& e2, const Event& e1) const {
            if (e1.x == e2.x) {
                const typename Polynomial<D>::CompareAt compare(e1.x);
                return compare(pieces[e1.id].polynomial, pieces[e2.id].polynomial);
            }
            return e1.x < e2.x;
        }
    };

    std::vector<PieceID> state;

    CompareEvent compare_event(pieces);
    std::priority_queue<Event, std::vector<Event>, CompareEvent> events(compare_event);

    for (PieceID id = 0; id < pieces.size(); ++id) {
        const PolynomialPiece<D>& piece = pieces[id];
        events.push({piece.interval.min, EventType::OPEN, id});
    }

    std::vector<PolynomialPiece<D>> result_pieces;
    // X of previous event
    double prev_x = std::numeric_limits<double>::min();
    // X of the start of the current piece
    double start_x = std::numeric_limits<double>::min();
    PieceID prev_open_id = INVALID_PIECE_ID;

    while (!events.empty()) {
        const auto event = events.top();
        events.pop();

        const typename Polynomial<D>::CompareAt compare(event.x);
        auto compare_pieces = [pieces, compare](PieceID p1, PieceID p2) -> bool {
            if (p1 == p2) return false;
            if (compare(pieces[p1].polynomial, pieces[p2].polynomial)) return true;
            if (compare(pieces[p2].polynomial, pieces[p1].polynomial)) return false;
            return p1 < p2;
        };

        if (event.type == EventType::OPEN) {
            const PolynomialPiece<D>& piece = pieces[event.id];

            // Add CLOSE event
            events.push({piece.interval.max, EventType::CLOSE, event.id});

            // Add SWAP events
            for (const auto other_id : state) {
                const PolynomialPiece<D>& other = pieces[other_id];
                const auto diff = piece.polynomial - other.polynomial;
                const auto roots = find_roots(diff);
                for (double root : roots) {
                    // Check that intersection is is contained in the pieces' interval.
                    // Exclude intersections at piece endpoints, as those are handled by OPEN / CLOSE events.
                    if (!piece.interval.contains_excl(root) || !other.interval.contains_excl(root)) {
                        continue;
                    }

                    // Check that the intersection actually swaps the order (and not just touches).
                    if (!diff.changes_sign_at(root)) {
                        continue;
                    }

                    // ID of the piece that is lowest after the event
                    PieceID lower_id = typename Polynomial<D>::CompareAt(root)(piece.polynomial, other.polynomial)
                                       ? event.id : other_id;

                    events.push({root, EventType::SWAP, lower_id});

//                    std::cout << "add swap " << event.id << ' ' << other_id << ' ' << lower_id << ' ' << root << '\n';
                }
            }

            // Insert piece into state at the correct position
            const auto geq_it = std::lower_bound(state.begin(), state.end(), event.id, compare_pieces);
            state.insert(geq_it, event.id);
        } else { // event.type == CLOSE || event.type == SWAP
            const typename Polynomial<D>::CompareAt compare_prev(prev_x);
            auto compare_pieces_prev = [pieces, compare_prev](PieceID p1, PieceID p2) -> bool {
                if (p1 == p2) return false;
                if (compare_prev(pieces[p1].polynomial, pieces[p2].polynomial)) return true;
                if (compare_prev(pieces[p2].polynomial, pieces[p1].polynomial)) return false;
                return p1 < p2;
            };

            // Iterator to the position of event.id in the state
            const auto it = std::lower_bound(state.begin(), state.end(), event.id, compare_pieces_prev);
            assert(*it == event.id);

            if (event.type == EventType::CLOSE) {
                // Remove piece from the state
                state.erase(it);
            } else { // event.type == SWAP
                assert(it != state.begin());
                std::iter_swap(it - 1, it);
            }
        }

        // Update result
        PieceID open_id = state.empty() ? INVALID_PIECE_ID : state.front();
        if (prev_open_id != open_id) {
            // Close previous piece
            if (
                // If we have already opened a piece
                prev_open_id != INVALID_PIECE_ID
                // And this would not result in an empty interval
                && start_x != event.x
            ) {
                result_pieces.emplace_back(Interval{start_x, event.x}, pieces[prev_open_id].polynomial);
            }

            // Open new piece
            if (open_id != INVALID_PIECE_ID) {
                start_x = event.x;
            }
        }

        prev_x = event.x;
        prev_open_id = open_id;

        // Verify state
        assert(std::is_sorted(state.begin(), state.end(), compare_pieces));

//        std::cout << "EVENT: x=" << event.x << " type=" << (size_t) event.type << "\n";
//        for (const auto id : state) {
//            std::cout << id << ' ' << pieces[id].polynomial(event.x) << ' ' << pieces[id] << '\n';
//        }
    }

    return {result_pieces};
}

template<size_t D>
void fast_lower_envelope(PiecewisePolynomial<D>& left, const PiecewisePolynomial<D>& right) {
    // Early return for trivial cases
    if (left.empty()) {
        left.pieces.assign(right.pieces.begin(), right.pieces.end());
        return;
    }
    if (right.empty()) {
        return;
    }
    if (left.interval().max == right.interval().min) {
        left.pieces.insert(left.pieces.end(), right.pieces.begin(), right.pieces.end());
        return;
    }

    // Verify that `left` is actually to the left of `right`
    assert(left.interval().min <= right.interval().min);
    assert(left.interval().max <= right.interval().max);
    // Verify that `left` and `right` have some overlap
    assert(left.interval().max > right.interval().min);

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

            // Shrink last piece of `left`
            left.pieces.back() = PolynomialPiece<D>({ left_back.interval.min, root }, left_back.polynomial);

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

    throw std::logic_error("unresolved case");
}

template<size_t dimension, Norm image_norm, Norm param_norm>
class CDTW
{
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

    // Dynamic program
    std::vector<std::vector<Entry>> in_functions;

    // Result
    PiecewisePolynomial<D> out_top;
    PiecewisePolynomial<D> out_right;

    // Internal methods
    PiecewisePolynomial<D> base_bottom(const Cell& cell) const;

public:
    CDTW(const Curve& curve1, const Curve& curve2);
};

template<size_t dimension, Norm image_norm, Norm param_norm>
CDTW<dimension, image_norm, param_norm>::CDTW(const Curve& curve1, const Curve& curve2) :
    in_functions(curve1.size() - 1, std::vector<Entry>(curve2.size() - 1))
{
    // TODO: Simplify cell structs?

    // Initialize bottom base in_functions
    {
        double c = 0;
        for (size_t i = 0; i + 1 < curve1.size(); ++i) {
            Cell cell(curve1[i], curve2[0], curve1[i + 1], curve2[1]);
            PiecewisePolynomial<D> f = base_bottom(cell) + c;
            in_functions[i][0].bottom = f;
            c = f.pieces.back().polynomial(f.pieces.back().interval.max); // TODO: Create and use back_value() instead

            std::cout << f << '\n';
        }
    }

    for (size_t i = 0; i + 1 < curve1.size(); ++i) {
        for (size_t j = 0; j + 1 < curve2.size(); ++j) {
            // In-functions along bottom/left are given
            // Compute out-functions along top/right
            // Depending on current row/column:
            //  - set out-functions as in-function for next cell
            //  - or append to out_top/out_right
        }
    }
}

//
// 1D + L1 image norm + L1 param norm
// TODO: Move to 1d-l1-l1.h
//

template<>
PiecewisePolynomial<2> CDTW<1, Norm::L1, Norm::L1>::base_bottom(const Cell& cell) const {
    // TODO: Clean this up; sx, sy, tx, ty are assumed to be coordinates in a space where (0,0) is the ellipse center
    double sx = cell.s.x - cell.mid.x;
    double sy = cell.s.y - cell.mid.y;
    double tx = cell.t.x - cell.mid.x;

    if (sx >= sy) {
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({-(sx * sx) / 2 + sx * sy, -sy, 1. / 2})}
        });
    } else if (sy >= tx) {
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})}
        });
    } else { // sx < sy < tx
        return PiecewisePolynomial<2>({
            {{sx + cell.mid.x, sy + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy, sy, -1. / 2})},
            {{sy + cell.mid.x, tx + cell.mid.x}, Polynomial<2>({(sx * sx) / 2 - sx * sy + sy * sy, -sy, 1. / 2})}
        });
    }
}
