#pragma once

#include <assert.h>

typedef double distance_t;

typedef arma::Row<distance_t>::fixed<2> Point;
typedef arma::Mat<distance_t> Points;

struct Line {
    [[deprecated]]
    distance_t slope;

    /**
     * y-intercept, or x-intercept if line is vertical
     */
    [[deprecated]]
    distance_t intercept;



    Line() = default;

    /**
     * Create line from two points that lie on the line.
     */
    Line(const Point& a, const Point& b) {
        assert(!arma::approx_equal(a, b, "absdiff", ABS_TOL));

        if (approx_equal(a(0), b(0))) { // Vertical line
            slope = a(1) < b(1) ? INFINITY : -INFINITY;
            intercept = a(0);
        } else {
            slope = (b(1) - a(1)) / (b(0) - a(0));
            intercept = a(1) - slope * a(0);
        }
    }

    /**
     * Create line from a point that lies on the line and its slope.
     */
    Line(const Point& p, distance_t slope) : slope(slope) {
        if (isVertical()) {
            intercept = p(0);
        } else {
            intercept = -slope * p(0) + p(1);
        }
    }

    inline bool isVertical() const {
        return std::isinf(slope);
    }

    inline bool isHorizontal() const {
        return slope == 0;
    }

    /**
     * Get Y coordinate for given X coordinate.
     */
    inline distance_t getY(distance_t x) const {
        if (isVertical()) {
            return NAN;
        }
        if (isHorizontal()) {
            return intercept;
        }
        return x * slope + intercept;
    }

    /**
     * Get X coordinate for given Y coordinate.
     */
    inline distance_t getX(distance_t y) const {
        if (isVertical()) {
            return intercept;
        }
        if (isHorizontal()) {
            return NAN;
        }
        return (y - intercept) / slope;
    }

    /**
     * Test whether the given point lies on this line.
     */
    bool includesPoint(const Point& point) const {
        if (isVertical()) {
            return approx_equal(intercept, point(0));
        } else {
            return approx_equal(getY(point(0)), point(1));
        }
    }

    /**
     * Find point at which the two given lines intersect.
     * ASSUMPTION: line1 and line2 are not parallel
     */
    friend Point intersect(const Line& line1, const Line& line2) {
        if (line2.isVertical()) {
            return intersect(line2, line1);
        }

        distance_t x;
        if (line1.isVertical()) {
            x = line1.intercept;
        } else {
            x = (line2.intercept - line1.intercept) / (line1.slope - line2.slope);
        }

        return {x, line2.slope * x + line2.intercept};
    }

    /**
     * Test whether the two given lines are perpendicular to each other.
     */
    friend bool isPerpendicular(const Line& line1, const Line& line2) {
        if (line2.isVertical()) {
            return isPerpendicular(line2, line1);
        }
        if (line1.isVertical()) {
            return line2.isHorizontal();
        }
        return approx_equal(line1.slope * line2.slope, -1.0);
    }
};
