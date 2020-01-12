#ifndef CODE_UTIL_H
#define CODE_UTIL_H

#include <cmath>

const double ABS_TOL = 1e-9;

template<class T>
bool approx_equal(T a, T b, double tol = ABS_TOL) {
    if (a == b) {
        return true;
    }

    if (std::isinf(a) || std::isinf(b)) { // && a != b
        return false;
    }

    return std::abs(a - b) <= tol;
}

template<class T>
bool approx_zero(T a, double tol = ABS_TOL) {
    return approx_equal(std::abs(a), 0.0, tol);
}

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

#endif //CODE_UTIL_H
