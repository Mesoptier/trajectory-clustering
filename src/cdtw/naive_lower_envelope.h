#ifndef TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H
#define TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H

#include <queue>
#include "PiecewisePolynomial.h"

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
//    std::cout << std::endl;

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
        CLOSE,
        OPEN,
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
                if (e1.type == e2.type) {
                    const typename Polynomial<D>::CompareAt compare(e1.x);
                    return compare(pieces[e1.id].polynomial, pieces[e2.id].polynomial);
                }
                return e1.type < e2.type;
            }
            return e1.x < e2.x;
        }
    };
    struct ComparePieceID
    {
        const std::vector<PolynomialPiece<D>>& pieces;
        const typename Polynomial<D>::CompareAt compare;
        ComparePieceID(const std::vector<PolynomialPiece<D>>& pieces, double x) : pieces(pieces), compare(x) {}

        bool operator()(PieceID p1, PieceID p2) const {
            if (p1 == p2) {
                // Exact same piece
                return false;
            }
            if (compare(pieces[p1].polynomial, pieces[p2].polynomial)) {
                // P1 is strictly lower
                return true;
            }
            if (compare(pieces[p2].polynomial, pieces[p1].polynomial)) {
                // P2 is strictly lower
                return false;
            }

            // Pieces have identical polynomials

            if (pieces[p1].interval.max != pieces[p2].interval.max) {
                // Prefer piece that ends later
                return pieces[p1].interval.max > pieces[p2].interval.max;
            }
            if (pieces[p1].interval.min != pieces[p2].interval.min) {
                // Prefer piece that starts earlier
                return pieces[p1].interval.min < pieces[p2].interval.min;
            }

            // Pieces have identical intervals and thus are entirely identical

            return p1 < p2;
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
    double prev_x = -std::numeric_limits<double>::infinity();
    // X of the start of the current piece
    double start_x = -std::numeric_limits<double>::infinity();
    PieceID prev_open_id = INVALID_PIECE_ID;

    while (!events.empty()) {
        const auto event = events.top();
        events.pop();

        ComparePieceID compare_pieces(pieces, event.x);

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
            ComparePieceID compare_pieces_prev(pieces, prev_x);

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
            if (prev_open_id != INVALID_PIECE_ID) {
                if (!approx_equal(start_x, event.x)) {
                    result_pieces.emplace_back(Interval{start_x, event.x}, pieces[prev_open_id].polynomial);
                }
            }

            // Open new piece
            if (open_id != INVALID_PIECE_ID) {
                if (!approx_equal(start_x, event.x)) {
                    start_x = event.x;
                }
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

#endif //TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H
