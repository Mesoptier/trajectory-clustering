#include "fast_lower_envelope.h"

template <size_t D>
PiecewisePolynomial<D> fast_lower_envelope_v2(const std::vector<PolynomialPiece<D>>& pieces) {
    if (pieces.size() == 1) 
        return PiecewisePolynomial<D>(pieces);

    auto left = fast_lower_envelope_v2(
        std::vector<PolynomialPiece<D>>(pieces.begin(), pieces.begin() + pieces.size() / 2)).pieces;
    auto right = fast_lower_envelope_v2(
        std::vector<PolynomialPiece<D>>(pieces.begin() + pieces.size() / 2, pieces.end())).pieces;

    return PiecewisePolynomial<D>(left);
}