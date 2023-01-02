#ifndef TRAJECTORY_CLUSTERING_INTERVAL_H
#define TRAJECTORY_CLUSTERING_INTERVAL_H

#include <iostream>
#include "../utils/util.h"

struct Interval_c {
    double min;
    double max;

    bool contains(double x) const;

    /**
     * Test whether x is contained in this interval (exclusive).
     * @param x
     * @return
     */
    bool contains_excl(double x) const;

    friend std::ostream& operator<<(std::ostream& os, const Interval_c& interval);

    bool operator==(const Interval_c& rhs) const;

    double interpolate(double t) const;

    Interval_c intersect(Interval_c& other);

    bool intersects(Interval_c& other);
};

template<>
bool approx_equal(const Interval_c& a, const Interval_c& b, double tol);

#endif //TRAJECTORY_CLUSTERING_INTERVAL_H
