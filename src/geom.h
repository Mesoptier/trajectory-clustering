#pragma once

#include <armadillo>
#include <assert.h>
#include "util.h"

using distance_t = double;

enum class Norm {
    L1,
    L2,
    LInf,
};

//
// Point
//

struct Point {
    distance_t x;
    distance_t y;

    Point() = default;
    Point(distance_t x, distance_t y) : x(x), y(y) {}

    Point& operator-=(const Point& point);
    Point operator-(const Point& point) const;
    Point& operator+=(const Point& point);
    Point operator+(const Point& point) const;
    Point& operator*=(distance_t mult);
    Point operator*(distance_t mult) const;
    Point& operator/=(distance_t distance);
    Point operator/(distance_t distance) const;

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;

    distance_t dist_sqr(const Point& point) const;
    distance_t dist(const Point& point) const;

    friend std::ostream& operator<<(std::ostream& out, const Point& p);
};

bool approx_equal(const Point& a, const Point& b);

/**
 * Computes perp dot product between two vectors.
 * See: http://mathworld.wolfram.com/PerpDotProduct.html
 */
distance_t perp(const Point& a, const Point& b);
distance_t dot(const Point& a, const Point& b);

distance_t norm(const Point& point, Norm p = Norm::L2);
Point normalise(const Point& point, Norm p = Norm::L2);

using Points = std::vector<Point>;

std::ostream& operator<<(std::ostream& out, const Points& points);

//
// Directions
//

// short for: backward-forward direction
enum class BFDirection {
    Backward = 0,
    Forward = 1,
};
BFDirection getMonotoneDirection(const Point& a, const Point &b);

struct MonotoneComparator {
    BFDirection direction;
    explicit MonotoneComparator(BFDirection direction): direction(direction) {}
    bool operator()(const Point& a, const Point& b);
};

struct Line
{
    Point origin;
    Point direction;

    Line() = default;

    Line(const Point& origin, const Point& direction)
        : origin(origin), direction(normalise(direction)) {}

    /**
     * Create line from two points that lie on the line.
     */
    static Line fromTwoPoints(const Point& a, const Point& b) {
        assert(!approx_equal(a, b));
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
        return approx_zero(direction.x);
    }

    inline bool isHorizontal() const {
        return approx_zero(direction.y);
    }

    /**
     * Get Y coordinate for given X coordinate.
     */
    inline distance_t getY(distance_t x) const {
        if (isVertical()) {
            return NAN;
        }
        return (x - origin.x) * direction.y / direction.x + origin.y;
    }

    /**
     * Get X coordinate for given Y coordinate.
     */
    inline distance_t getX(distance_t y) const {
        if (isHorizontal()) {
            return NAN;
        }
        return (y - origin.y) * direction.x / direction.y + origin.x;
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
        const Point perp1 = {-line1.direction.y, line1.direction.x};
        const auto t2 =
            dot(perp1, line1.origin - line2.origin) / dot(perp1, line2.direction);
        return line2(t2);
    }

    friend bool isParallel(const Line& line1, const Line& line2) {
        return approx_equal(std::abs(dot(line1.direction, line2.direction)), 1.0);
    }

    friend bool isSameDirection(const Line& line1, const Line& line2) {
        return approx_equal(dot(line1.direction, line2.direction), 1.0);
    }

    friend bool isOppositeDirection(const Line& line1, const Line& line2) {
        return approx_equal(dot(line1.direction, line2.direction), -1.0);
    }

    /**
     * Test whether the two given lines are perpendicular to each other.
     */
    friend bool isPerpendicular(const Line& line1, const Line& line2) {
        return approx_zero(dot(line1.direction, line2.direction));
    }
};
