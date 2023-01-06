#include "PiecewisePolynomial.h"

template <size_t D>
std::vector<PolynomialPiece<D>> reduce_list(std::vector<PolynomialPiece<D>> pieces) {
    auto output = std::vector<PolynomialPiece<D>>();

    for (int i = 0; i < pieces.size(); ++i) {
        auto piece = pieces[i];
        
        // if (i > 0)
        //     assert(approx_equal(pieces[i-1].interval.max, piece.interval.min));

        if (!output.empty() && output.back().polynomial == piece.polynomial) {
            output.back().interval.max = piece.interval.max;
        } else {
            output.push_back(piece);
        }
    }
    return output;
}

template <size_t D>
bool are_equal(std::vector<PolynomialPiece<D>> set1, std::vector<PolynomialPiece<D>> set2) {
    if (set1.size() != set2.size())
        return false;

    for (int i = 0; i < set1.size(); ++i) {
        auto p1 = set1[i];
        auto p2 = set2[i];

        if (!approx_equal(p1.interval.min, p2.interval.min) || 
        !approx_equal(p1.interval.max, p2.interval.max))
            return false;

        if (p1.polynomial != p2.polynomial) {
            std::cout << "polynomials not equal\n";
            return false;
        }
    }

    return true;
}

template <size_t D>
std::vector<PolynomialPiece<2>> read_pieces(std::string filename) {
    std::ifstream file(filename);
    std::string line;

    std::vector<PolynomialPiece<2>> output = std::vector<PolynomialPiece<2>>();

    distance_t a, b, c, min, max;
    while (file >> a >> b >> c >> min >> max && std::getline(file, line)) {
        output.push_back(PolynomialPiece<2>({min, max},
        Polynomial<2>({{a, b, c}})));
    }

    return output;
}

/**
 * @brief Computes the lower envelope of the the given polynomial pieces
 * using a divide and conquer algorithm. The output is in left-to-right order
 * with respect to the intervals of the polynomials pieces.
 * 
 * @tparam D degree of polynomals
 * @param pieces input polynomials
 * @return PiecewisePolynomial<D> lower envelope in left-to-right order
 */
template <size_t D>
PiecewisePolynomial<D> fast_lower_envelope_v2(const std::vector<PolynomialPiece<D>>& pieces) {
    if (pieces.empty())
        return {};
    
    if (pieces.size() == 1) 
        return PiecewisePolynomial<D>(pieces);

    auto output_pieces = std::vector<PolynomialPiece<D>>();

    // Recursively compute the lower envelopes of the left and right halves.
    auto left = fast_lower_envelope_v2(
        std::vector<PolynomialPiece<D>>(pieces.begin(), pieces.begin() + pieces.size() / 2)
    ).pieces;
    auto right = fast_lower_envelope_v2(
        std::vector<PolynomialPiece<D>>(pieces.begin() + pieces.size() / 2, pieces.end())
    ).pieces;

    // Reverse the order of left and right so the left-most piece
    // can be added/removed in constant time
    std::reverse(left.begin(), left.end());
    std::reverse(right.begin(), right.end());

    while (!left.empty() || !right.empty()) {
        if (left.empty()) {
            // When left is empty, add the remaing pieces from right
            // to the output list in left-to-right order
            for (int i = right.size() - 1; i >= 0; --i) {
                auto p = right[i];
                output_pieces.push_back(p);
            }
            right = std::vector<PolynomialPiece<D>>();
        } else if (right.empty()) {
            // similarly when right is empty...
            for (int i = left.size() - 1; i >= 0; --i) {
                auto p = left[i];
                output_pieces.push_back(p);
            }
            left = std::vector<PolynomialPiece<D>>();
        } else {
            // When left and right are non-empty, remove the left-most
            // piece from each list
            auto left_piece = left.back();
            auto right_piece = right.back();
            left.pop_back();
            right.pop_back();

            std::vector<PolynomialPiece<D>> new_pieces = std::vector<PolynomialPiece<D>>();
            // break_points will store the intersection points
            // of left_piece and right_piece as well as the endpoitns of the
            // intersection of their intervals.
            std::vector<double> break_points;

            if (left_piece.interval.intersects(right_piece.interval)) {

                if (left_piece.interval.min < right_piece.interval.min &&
                !approx_equal(left_piece.interval.min, right_piece.interval.min)) {
                    // There is a subinterval of left_piece.interval which is 
                    // completely to the left of right_piece.interval
                    // A polynomial piece definied on this interval with polynomial
                    // left_piece.polynomial can be added to the output.
                    PolynomialPiece<D> new_piece = PolynomialPiece<D>(
                        {left_piece.interval.min, right_piece.interval.min}, left_piece.polynomial
                    );
                    new_piece.history = left_piece.history;
                    new_pieces.push_back(new_piece);
                    
                } else if(right_piece.interval.min < left_piece.interval.min &&
                !approx_equal(right_piece.interval.min, left_piece.interval.min)) {
                    // Analgously handle the case where right_piece.interval is left-most...
                    PolynomialPiece<D> new_piece = PolynomialPiece<D>(
                        {right_piece.interval.min, left_piece.interval.min}, right_piece.polynomial
                    );
                    new_piece.history = right_piece.history;
                    new_pieces.push_back(new_piece);
                }

                auto intersections = find_intersections(left_piece.polynomial, right_piece.polynomial);
                Interval_c interval = left_piece.interval.intersect(right_piece.interval);
                break_points = {interval.min, interval.max};
                for (double p: intersections) {
                    // Only necessary to consider intersection points which lie inside
                    // interval since outside this interval left_piece and right_piece
                    // cannot both be valid.
                    if (interval.contains_excl(p) &&
                    !approx_equal(p, interval.min) && !approx_equal(p, interval.max))
                        break_points.push_back(p);
                }
                std::sort(break_points.begin(), break_points.end());

                    for (int i = 1; i < break_points.size(); ++i) {
                        double mid = (break_points[i-1] + break_points[i]) / 2;
                        Polynomial<D> new_poly;
                        // between each pair of consecutive break_points, 
                        // set new_poly to be the polynomial which attains 
                        // a smaller value on this interval.
                        // It is sufficient compare both polynomials at mid
                        // since there can be no intersections inside the interval
                        // {break_points[i-1], break_points[i]}
                        if (left_piece.polynomial(mid) < right_piece.polynomial(mid)) {
                            new_poly = left_piece.polynomial;
                        } else {
                            new_poly = right_piece.polynomial;
                        }

                        auto new_piece = PolynomialPiece<D>(
                                    {break_points[i-1], break_points[i]},
                                    new_poly
                                );
                        new_piece.history = left_piece.polynomial(mid) < right_piece.polynomial(mid) ?
                        left_piece.history : right_piece.history;
                    
                        // We exculde polynomials defined on a single point.
                        if (output_pieces.size() == 0 ||
                                output_pieces.back().interval.max < break_points[i])
                            new_pieces.push_back(new_piece);
                    }

                // The lower envelope has now been computed on all points on which both
                // left_piece and right_piece are defined.
                // There may be an interval to the right of this common interval where only
                // one of the two polynomials are defined. This piece needs to be return to
                // to its originating vector so that it can be compared with other pieces from
                // the other vector.
                if (left_piece.interval.max < right_piece.interval.max && 
                !approx_equal(left_piece.interval.max, right_piece.interval.max)) {
                    // return the remaineder of right_piece to its vector
                    right.push_back(
                        PolynomialPiece<D>(
                            {left_piece.interval.max, right_piece.interval.max},
                            right_piece.polynomial
                        )
                    );
                    right.back().history = right_piece.history;
                } else if (right_piece.interval.max < left_piece.interval.max && 
                !approx_equal(right_piece.interval.max, left_piece.interval.max)) {
                    // return the remaineder of left_piece to its vector
                    left.push_back(
                        PolynomialPiece<D>(
                            {right_piece.interval.max, left_piece.interval.max},
                            left_piece.polynomial
                        )
                    );
                    left.back().history = left_piece.history;
                }

            } else {
                // If the intervals do not intersect at all, the piece with
                // the left-most interval can be appended to the output
                // while the other piece must be returned to its vector
                if (left_piece.interval.min < right_piece.interval.min) {
                    new_pieces.push_back(left_piece);
                    right.push_back(right_piece);
                } else {
                    new_pieces.push_back(right_piece);
                    left.push_back(left_piece);
                }
            }

            new_pieces = reduce_list(new_pieces);
            
            for (auto p: new_pieces) {
                if (!approx_equal(p.interval.min, p.interval.max))
                    output_pieces.push_back(p);
            }
            output_pieces = reduce_list(output_pieces);
        }
    }

    if (output_pieces.size() == 0)
        return {};

    output_pieces = reduce_list(output_pieces);
    return PiecewisePolynomial<D>(output_pieces);
}