#ifndef TRAJECTORY_CLUSTERING_INTERVAL_H
#define TRAJECTORY_CLUSTERING_INTERVAL_H

#include <iostream>

struct Interval {
    double min;
    double max;

    bool contains(double x) const {
        return min <= x && x <= max;
    }

    /**
     * Test whether x is contained in this interval (exclusive).
     * @param x
     * @return
     */
    bool contains_excl(double x) const {
        return min < x && x < max;
    }

    friend std::ostream& operator<<(std::ostream& os, const Interval& interval) {
        os << interval.min << " <= x <= " << interval.max;
        return os;
    }
};

#endif //TRAJECTORY_CLUSTERING_INTERVAL_H
