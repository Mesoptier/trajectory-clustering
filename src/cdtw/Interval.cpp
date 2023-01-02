#include "Interval.h"
#include "../utils/util.h"

template<>
bool approx_equal(const Interval_c& a, const Interval_c& b, double tol) {
    return approx_equal(a.min, b.min, tol) && approx_equal(a.max, b.max, tol);
}

bool Interval_c::contains(double x) const {
        return min <= x && x <= max ||
        approx_equal(x, max) || approx_equal(x, min);
}

bool Interval_c::contains_excl(double x) const {
    return min + ABS_TOL <= x && x <= max - ABS_TOL;
}

std::ostream& operator<<(std::ostream& os, const Interval_c& interval) {
    os << interval.min << " <= x <= " << interval.max;
    return os;
}

bool Interval_c::operator==(const Interval_c& rhs) const {
    return min == rhs.min && max == rhs.max;
}

double Interval_c::interpolate(double t) const {
    return min * (1 - t) + max * t;
}

Interval_c Interval_c::intersect(Interval_c& other) {
    
    double left = std::max(min, other.min);
    double right = std::min(max, other.max);

    return { left, right };
}

bool Interval_c::intersects(Interval_c& other) {
    return (contains_excl(other.min) || contains_excl(other.max))
    || (other.contains_excl(min) || other.contains_excl(max)) ||
    (approx_equal(min, other.min) && approx_equal(max, other.max));
}