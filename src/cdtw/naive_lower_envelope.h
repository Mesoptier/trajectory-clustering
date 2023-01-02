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
    if (pieces.empty()) {
        return {};
    }
    auto start_time = std::chrono::high_resolution_clock::now();
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
        SWAP,
        OPEN,
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
            if (e1.x != e2.x) return e1.x < e2.x;
            if (!(e1.type == e2.type)) return e1.type < e2.type;
            if (e1.id != e2.id) return e1.id < e2.id;
            const typename Polynomial<D>::CompareAt compare(e1.x);
            return compare(pieces[e1.id].polynomial, pieces[e2.id].polynomial);
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
        // double x = event.x;
        // double height = pieces[event.id](x);
        
        // std::cout << event.x << ", " <<  pieces[event.id] << std::endl;
        
        ComparePieceID compare_pieces(pieces, event.x);

        if (event.type == EventType::OPEN) {
            const PolynomialPiece<D>& piece = pieces[event.id];

            // Add CLOSE event
            events.push({piece.interval.max, EventType::CLOSE, event.id});

            // Add SWAP events
            for (const auto other_id : state) {
                const PolynomialPiece<D>& other = pieces[other_id];
                double overlap_min = std::max(piece.interval.min, other.interval.min);
                double overlap_max = std::min(piece.interval.max, other.interval.max);
                if (overlap_max < overlap_min) {
                    continue;
                }

                const auto diff = piece.polynomial - other.polynomial;
                const auto roots = find_roots(diff);
                auto debug_diff = diff.to_string();
                for (double root : roots) {
                    // Check that intersection is is contained in the pieces' interval.
                    // Exclude intersections at piece endpoints, if those would be handled by OPEN / CLOSE events.
                    if (overlap_min <= root && root <= overlap_min + ABS_TOL) {
                        if (approx_zero(diff(overlap_min))) {
                            continue;
                        }
                    }
                    if (overlap_max - ABS_TOL <= root && root <= overlap_max) {
                        if (approx_zero(diff(overlap_max))) {
                            continue;
                        }
                    }
                    if (root < overlap_min || overlap_max < root) {
                        continue;
                    }

                    // std::cout << diff(root) << std::endl;
                    // assert(approx_zero(diff(root)));

                    // Check that the intersection actually swaps the order (and not just touches).
                    if (!diff.changes_sign_at(root)) {
                        continue;
                    }

                    // ID of the piece that is lowest after the event
                    PieceID lower_id = typename Polynomial<D>::CompareAt(root)(piece.polynomial, other.polynomial)
                                       ? event.id : other_id;

                    events.push({root, EventType::SWAP, lower_id});
                }
            }

            // Insert piece into state at the correct position
            const auto geq_it = std::lower_bound(state.begin(), state.end(), event.id, compare_pieces);
            state.insert(geq_it, event.id);
        } else { // event.type == CLOSE || event.type == SWAP
            ComparePieceID compare_pieces_prev(pieces, prev_x);

            // Iterator to the position of event.id in the state
//            auto it = std::lower_bound(state.begin(), state.end(), event.id, compare_pieces_prev);
            auto it = std::find_if(state.begin(), state.end(), [event](auto id) { return id == event.id; });
            assert(*it == event.id);

            if (event.type == EventType::CLOSE) {
                // Remove piece from the state
                state.erase(it);
            } else { // event.type == SWAP
//                assert(it != state.begin());
//                std::iter_swap(it - 1, it);
//                --it;
//
//                // Resolve SWAP events in the same point
//                while (events.top().type == EventType::SWAP && events.top().x == event.x) {
//                    const auto next_event = events.top();
//                    events.pop();
//
//                    // The piece with next_event.id is either at `it` or above it in the state
//                    while (*it != next_event.id) {
//                        ++it;
//                        assert(it != state.end());
//                    }
//
//                    assert(it != state.begin());
//                    std::iter_swap(it - 1, it);
//                    --it;
//                }
            }
        }

        // TODO: Figure out how to use the SWAP events again
        std::sort(state.begin(), state.end(), compare_pieces);

        // Update result
        PieceID open_id = state.empty() ? INVALID_PIECE_ID : state.front();
        if (prev_open_id != open_id) {
            // Close previous piece
            if (prev_open_id != INVALID_PIECE_ID) {
                if (!approx_equal(start_x, event.x)) {
                    if (!result_pieces.empty() && approx_equal(result_pieces.back().polynomial, pieces[prev_open_id].polynomial)
                    && result_pieces.back().history.path_type == pieces[prev_open_id].history.path_type) {
                        if (approx_equal(result_pieces.back().test_value, 0.67677669529663675))
                            std::cout << "hi\n";
                        result_pieces.back().interval.max = event.x;
                    } else {
                        result_pieces.emplace_back(Interval_c{start_x, event.x}, pieces[prev_open_id].polynomial);
                        result_pieces.back().history = pieces[prev_open_id].history;
                        result_pieces.back().test_value = result_pieces.back().polynomial(result_pieces.back().interval.max);
                        if (approx_equal(result_pieces.back().test_value, 0.67677669529663675))
                            std::cout << "hi\n";
                    }
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
    }

    // for (int i = 1; i < result_pieces.size(); ++i) {
    //         auto l = result_pieces[i-1];
    //         auto r = result_pieces[i];
    //         assert(approx_equal(l.interval.max, r.interval.min));
    //         auto l_max = l.polynomial(l.interval.max);
    //         auto r_min = r.polynomial(r.interval.min);
    //         if (!(approx_equal(l_max, r_min)))
    //             assert(approx_equal(l_max, r_min));
    //     }



    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    return {result_pieces};
}

#endif //TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H
