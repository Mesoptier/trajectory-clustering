#ifndef TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H
#define TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H


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

#endif //TRAJECTORY_CLUSTERING_NAIVE_LOWER_ENVELOPE_H
