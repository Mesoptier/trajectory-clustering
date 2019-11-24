#ifndef CODE_GEOM_H
#define CODE_GEOM_H

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
        // ASSUMPTION: a != b

        if (approx_equal(a(0), b(0))) { // Vertical line
            slope = a(1) < b(1) ? INFINITY : -INFINITY;
            intercept = a(0);
        } else {
            slope = (b(1) - a(1)) / (b(0) - a(0));
            intercept = a(1) - slope * a(0);
        }
    }

    bool isVertical() const {
        return std::isinf(slope);
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
