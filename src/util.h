#ifndef CODE_UTIL_H
#define CODE_UTIL_H

#include <cmath>

const double ABS_TOL = 1e-13;

template<class T>
bool approx_equal(T a, T b, double tol = ABS_TOL) {
    if (std::signbit(a) != std::signbit(b)) {
        return false;
    }

    if (std::isinf(a) && a == b) {
        return true;
    }

    if (std::isinf(a) != std::isinf(b)) {
        return false;
    }

    return std::abs(a - b) <= tol;
}

template<class T>
bool approx_zero(T a, double tol = ABS_TOL) {
    return approx_equal(std::abs(a), 0.0, tol);
}

#endif //CODE_UTIL_H
