#include "Interval.h"

template<>
bool approx_equal(const Interval& a, const Interval& b, double tol) {
    return approx_equal(a.min, b.min, tol) && approx_equal(a.max, b.max, tol);
}