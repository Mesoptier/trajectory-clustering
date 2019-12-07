#pragma once

#include <armadillo>
#include <assert.h>
#include "util.h"

typedef double distance_t;

typedef arma::Row<distance_t>::fixed<2> Point;
typedef arma::Mat<distance_t> Points;
typedef std::vector<Point> PointsList;

/**
 * Computes perp dot product between two vectors.
 * See: http://mathworld.wolfram.com/PerpDotProduct.html
 */
distance_t perp(const Point& a, const Point& b);

struct MonotoneComparator {
    enum Direction {
        LowerFirst,
        HigherFirst,
    };

    Direction direction;

    explicit MonotoneComparator(Direction direction): direction(direction) {}

    bool operator()(const Point& a, const Point& b) {
        return (direction == LowerFirst)
               ? a(0) < b(0) + ABS_TOL && a(1) < b(1) + ABS_TOL
               : a(0) + ABS_TOL > b(0) && a(1) + ABS_TOL > b(1);
    }

    static Direction getDirection(const Point& a, const Point& b) {
        #ifndef NDEBUG
        if (arma::approx_equal(a, b, "absdiff", ABS_TOL)) {
            throw std::logic_error("Points are equal, and therefore not monotone");
        }
        #endif

        if (MonotoneComparator(LowerFirst)(a, b)) {
            return LowerFirst;
        }
        if (MonotoneComparator(HigherFirst)(a, b)) {
            return HigherFirst;
        }

        std::stringstream error;
        error << "Points (" << a(0) << "," << a(1) << ") and (" << b(0) << "," << b(1) << ") are not monotone";
        throw std::logic_error(error.str());
    }
};

struct Line
{
    Point origin;
    Point direction;

    Line() = default;

    Line(const Point& origin, const Point& direction)
        : origin(origin), direction(arma::normalise(direction)) {}

    /**
     * Create line from two points that lie on the line.
     */
    static Line fromTwoPoints(const Point& a, const Point& b) {
        assert(!arma::approx_equal(a, b, "absdiff", ABS_TOL));
        return {a, b - a};
    }

    /**
     * Create line from a point that lies on the line and its slope.
     */
    static Line fromPointAndSlope(Point p, distance_t slope) {
        if (slope == INFINITY) {
            return {p, {0, 1}};
        }
        if (slope == -INFINITY) {
            return {p, {0, -1}};
        }
        return {p, {1, slope}};
    }

    Point operator()(distance_t t) const {
        return origin + direction * t;
    }

    inline bool isVertical() const {
        return approx_zero(direction(0));
    }

    inline bool isHorizontal() const {
        return approx_zero(direction(1));
    }

    /**
     * Get Y coordinate for given X coordinate.
     */
    inline distance_t getY(distance_t x) const {
        if (isVertical()) {
            return NAN;
        }
        return (x - origin(0)) * direction(1) / direction(0) + origin(1);
    }

    /**
     * Get X coordinate for given Y coordinate.
     */
    inline distance_t getX(distance_t y) const {
        if (isHorizontal()) {
            return NAN;
        }
        return (y - origin(1)) * direction(0) / direction(1) + origin(0);
    }

    /**
     * Test whether the given point lies on this line.
     */
    bool includesPoint(const Point& point) const {
        return approx_zero(perp(origin - point, direction));
    }

    /**
     * Find point at which the two given lines intersect.
     * ASSUMPTION: line1 and line2 are not parallel
     */
    friend Point intersect(const Line& line1, const Line& line2) {
        const Point perp1 = {-line1.direction(1), line1.direction(0)};
        const auto t2 =
            arma::dot(perp1, line1.origin - line2.origin) / arma::dot(perp1, line2.direction);
        return line2(t2);
    }

    friend bool isParallel(const Line& line1, const Line& line2) {
        return approx_equal(std::abs(arma::dot(line1.direction, line2.direction)), 1.0);
    }

    friend bool isSameDirection(const Line& line1, const Line& line2) {
        return approx_equal(arma::dot(line1.direction, line2.direction), 1.0);
    }

    friend bool isOppositeDirection(const Line& line1, const Line& line2) {
        return approx_equal(arma::dot(line1.direction, line2.direction), -1.0);
    }

    /**
     * Test whether the two given lines are perpendicular to each other.
     */
    friend bool isPerpendicular(const Line& line1, const Line& line2) {
        return approx_zero(arma::dot(line1.direction, line2.direction));
    }
};
