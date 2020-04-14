#include "gtest/gtest.h"
#include "../src/cdtw/naive_lower_envelope.h"

void AssertIntervalEqual(const Interval& actual, const Interval& expected) {
    ASSERT_DOUBLE_EQ(actual.min, expected.min);
    ASSERT_DOUBLE_EQ(actual.max, expected.max);
}

template<size_t D>
void AssertPolynomialEqual(const Polynomial<D>& actual, const Polynomial<D>& expected) {
    for (size_t d = 0; d <= D; ++d) {
        ASSERT_DOUBLE_EQ(actual.coefficients.at(d), expected.coefficients.at(d));
    }
}

template<size_t D>
void AssertPiecewisePolynomialEqual(const PiecewisePolynomial<D>& actual, const PiecewisePolynomial<D>& expected) {
    ASSERT_EQ(actual.pieces.size(), expected.pieces.size());
    for (size_t i = 0; i < actual.pieces.size(); ++i) {
        const PolynomialPiece<D>& actual_piece = actual.pieces.at(i);
        const PolynomialPiece<D>& expected_piece = expected.pieces.at(i);

        AssertIntervalEqual(actual_piece.interval, expected_piece.interval);
        AssertPolynomialEqual(actual_piece.polynomial, expected_piece.polynomial);
    }
}

// TODO, test the following:
//  - Pieces with identical polynomials
//  - Multiple pieces intersecting in a single point
//  - Pieces intersecting at their start-/endpoints
//  - Possible floating point issues
//  - More?

TEST(NaiveLowerEnvelopeTest, Basic) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({0, 0, 1})},
        {{0, 2}, Polynomial<2>({4, -4, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        {{0, 1}, Polynomial<2>({0, 0, 1})},
        {{1, 2}, Polynomial<2>({4, -4, 1})},
    });
    AssertPiecewisePolynomialEqual(actual, expected);
}