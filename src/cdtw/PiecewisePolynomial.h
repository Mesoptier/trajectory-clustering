#ifndef TRAJECTORY_CLUSTERING_PIECEWISEPOLYNOMIAL_H
#define TRAJECTORY_CLUSTERING_PIECEWISEPOLYNOMIAL_H

#include "Interval.h"
#include "Polynomial.h"

template<size_t D>
struct PolynomialPiece
{
    Interval interval;
    Polynomial<D> polynomial;

    PolynomialPiece(const Interval& interval, const Polynomial<D>& polynomial) :
        interval(interval), polynomial(polynomial) {}

    PolynomialPiece<D> translate(double cx) const {
        auto result = *this;
        result.interval.min += cx;
        result.interval.max += cx;
        result.polynomial = result.polynomial.translate_xy(cx, 0);
        return result;
    }

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
    // Equality operators
    //

    bool operator==(const PolynomialPiece& rhs) const {
        return interval == rhs.interval &&
            polynomial == rhs.polynomial;
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
            assert(approx_equal(p1.interval.max, p2.interval.min));
            // Verify that pieces connect
            assert(approx_equal(p1.polynomial(p1.interval.max), p2.polynomial(p2.interval.min)));
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

    PiecewisePolynomial<D> translate(double cx) const {
        auto result = *this;
        for (auto& piece : result.pieces) {
            piece = piece.translate(cx);
        }
        return result;
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
    // Equality operators
    //

    bool operator==(const PiecewisePolynomial& rhs) const {
        return pieces == rhs.pieces;
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

#endif //TRAJECTORY_CLUSTERING_PIECEWISEPOLYNOMIAL_H
