#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <functional>

double const ABS_TOL = 1e-9;

template<class T>
bool approx_equal(T const& a, T const& b, double tol = ABS_TOL) {
    if (a == b)
        return true;
    if (std::isinf(a) || std::isinf(b)) // && a != b
        return false;
    return std::abs(a - b) <= tol;
}

template<class T>
bool approx_zero(T a, double tol = ABS_TOL) {
    return approx_equal(std::abs(a), 0.0, tol);
}

template <class T>
inline void hash_combine(std::size_t& seed, T const& v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}
#endif //CODE_UTIL_H
