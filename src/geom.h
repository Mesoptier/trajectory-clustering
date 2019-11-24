#ifndef CODE_GEOM_H
#define CODE_GEOM_H

#include <assert.h>

typedef double distance_t;

typedef arma::Row<distance_t>::fixed<2> Point;

struct Line {
    distance_t slope;

    /**
     * y-intercept, or x-intercept if line is vertical
     */
    distance_t intercept;

    Line() = default;

    Line(Point a, Point b) {
        assert(!arma::approx_equal(a, b, "absdiff", ABS_TOL));

        if (approx_equal(a(0), b(0))) { // Vertical line
            slope = a(1) < b(1) ? INFINITY : -INFINITY;
            intercept = a(0);
        } else {
            slope = (b(1) - a(1)) / (b(0) - a(0));
            intercept = a(1) - slope * a(0);
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
     * Tests whether the given point lies on this line.
     */
    bool includesPoint(const Point& point) const {
        if (isVertical()) {
            return approx_equal(intercept, point(0));
        } else {
            return approx_equal(getY(point(0)), point(1));
        }
    }

    friend Point intersect(const Line& line1, const Line& line2) {
        // ASSUMPTION: line1 and line2 are not parallel

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
};

#endif //CODE_GEOM_H
