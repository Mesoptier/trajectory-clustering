#ifndef TRAJECTORY_CLUSTERING_INTERVAL_H
#define TRAJECTORY_CLUSTERING_INTERVAL_H

#include <iostream>
#include "../util.h"

struct Interval {
    double min;
    double max;

    bool contains(double x) const {
        return min <= x && x <= max ||
        approx_equal(x, max) || approx_equal(x, min);
    }

    /**
     * Test whether x is contained in this interval (exclusive).
     * @param x
     * @return
     */
    bool contains_excl(double x) const {
        return min + ABS_TOL <= x && x <= max - ABS_TOL;
    }

    friend std::ostream& operator<<(std::ostream& os, const Interval& interval) {
        os << interval.min << " <= x <= " << interval.max;
        return os;
    }

    bool operator==(const Interval& rhs) const {
        return min == rhs.min && max == rhs.max;
    }

    double interpolate(double t) const {
        return min * (1 - t) + max * t;
    }

    Interval intersect(Interval& other) {
        
        double left = std::max(min, other.min);
        double right = std::min(max, other.max);

        return { left, right };
    }
};

template<>
bool approx_equal(const Interval& a, const Interval& b, double tol);

#endif //TRAJECTORY_CLUSTERING_INTERVAL_H
