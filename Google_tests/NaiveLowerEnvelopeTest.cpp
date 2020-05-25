#include "gtest/gtest.h"
#include "../src/cdtw/naive_lower_envelope.h"

constexpr double EPS = 1e-11;

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
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, IdenticalPieces) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({3, 2, 1})},
        {{0, 2}, Polynomial<2>({3, 2, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        {{0, 2}, Polynomial<2>({3, 2, 1})},
    });
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, IdenticalPolynomials) {
    //
    // Different start-/endpoints
    //
    { // Order: left, right
        const std::vector<PolynomialPiece<2>> pieces {
            {{0, 2}, Polynomial<2>({3, 2, 1})},
            {{1, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 1}, Polynomial<2>({3, 2, 1})},
            {{1, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
    { // Order: right, left
        const std::vector<PolynomialPiece<2>> pieces {
            {{1, 3}, Polynomial<2>({3, 2, 1})},
            {{0, 2}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 1}, Polynomial<2>({3, 2, 1})},
            {{1, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }

    //
    // Same startpoints
    //
    { // Order: left, right
        const std::vector<PolynomialPiece<2>> pieces {
            {{0, 2}, Polynomial<2>({3, 2, 1})},
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
    { // Order: right, left
        const std::vector<PolynomialPiece<2>> pieces {
            {{0, 3}, Polynomial<2>({3, 2, 1})},
            {{0, 2}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
    { // Order: left, right
        const std::vector<PolynomialPiece<2>> pieces {
            {{0 + EPS, 2}, Polynomial<2>({3, 2, 1})},
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
    { // Order: left, right
        const std::vector<PolynomialPiece<2>> pieces {
            {{0, 2}, Polynomial<2>({3, 2, 1})},
            {{0 + EPS, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }

    //
    // Same endpoints
    //
    { // Order: left, right
        const std::vector<PolynomialPiece<2>> pieces {
            {{0, 3}, Polynomial<2>({3, 2, 1})},
            {{1, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
    { // Order: right, left
        const std::vector<PolynomialPiece<2>> pieces {
            {{1, 3}, Polynomial<2>({3, 2, 1})},
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        };
        const auto actual = naive_lower_envelope(pieces);
        const PiecewisePolynomial<2> expected({
            {{0, 3}, Polynomial<2>({3, 2, 1})},
        });
        ASSERT_EQ(actual, expected);
    }
}

TEST(NaiveLowerEnvelopeTest, SlightlyDifferentEndpoints) {
    // The first function is significantly higher than the second at all points, but its endpoint is slightly later and
    // so one might expect it to show up in the result. But this would cause disconnected pieces.
    // So we don't want to see the first function in the result.
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2 + EPS}, Polynomial<2>({5, 2, 1})},
        {{0, 2}, Polynomial<2>({3, 2, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        // TODO: We might want for the result to be in range {0, 2 + EPS}
        {{0, 2}, Polynomial<2>({3, 2, 1})},
    });
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, SlightlyDifferentStartpoints) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({5, 2, 1})},
        {{0 + EPS, 2}, Polynomial<2>({3, 2, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        {{0, 2}, Polynomial<2>({3, 2, 1})},
    });
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, IdenticalPolynomialsHigherAfterIntersection) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({1, -2, 1})},
        {{0, 2}, Polynomial<2>({1, -2, 1})},
        {{0, 2}, Polynomial<2>({1, -2, 1})},
        {{1, 3}, Polynomial<2>({4, -4, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        {{0, 1.5}, Polynomial<2>({1, -2, 1})},
        {{1.5, 3}, Polynomial<2>({4, -4, 1})},
    });
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, IdenticalPolynomialsLowerAfterIntersection) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({1, -2, 1})},
        {{1, 3}, Polynomial<2>({4, -4, 1})},
        {{1, 3}, Polynomial<2>({4, -4, 1})},
        {{1, 3}, Polynomial<2>({4, -4, 1})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        {{0, 1.5}, Polynomial<2>({1, -2, 1})},
        {{1.5, 3}, Polynomial<2>({4, -4, 1})},
    });
    ASSERT_EQ(actual, expected);
}

TEST(NaiveLowerEnvelopeTest, PieceEndsSlightlyBeforeIntersection) {
    const std::vector<PolynomialPiece<2>> pieces {
        {{0, 2}, Polynomial<2>({0, 1, 0})},
        {{0, 1 - EPS}, Polynomial<2>({0, 1, 0})},
        {{0, 2}, Polynomial<2>({2, -1, 0})},
    };
    const auto actual = naive_lower_envelope(pieces);
    const PiecewisePolynomial<2> expected({
        // TODO: Instead of `1-EPS` we want the breakpoint to be at `1` exactly
        {{0, 1-EPS}, Polynomial<2>({0, 1, 0})},
        {{1-EPS, 2}, Polynomial<2>({2, -1, 0})},
    });
    ASSERT_EQ(actual, expected);
}