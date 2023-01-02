#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <functional>

// double const ABS_TOL = 1e-12;
// double const ABS_TOL = 1e-9;
double const ABS_TOL = 1e-6;
// double const ABS_TOL = 1e-4;


// int lower_enelope_time;
// int total_cdtw_time;
// int solve_integral_time;


template<class T>
bool approx_equal(T const& a, T const& b, double tol = ABS_TOL) {
    if (a == b)
        return true;
    if (std::isinf(a) || std::isinf(b)) // && a != b
        return false;
    return std::abs(a - b) <= tol;
}

template<class T>
bool approx_equal(const std::vector<T>& a, const std::vector<T>& b, double tol = ABS_TOL) {
    if (a.size() != b.size()) {
        return false;
    }

    for (size_t i = 0; i < a.size(); ++i) {
        if (!approx_equal(a[i], b[i], tol)) {
            return false;
        }
    }

    return true;
}

template<class T, size_t N>
bool approx_equal(const std::array<T, N>& a, const std::array<T, N>& b, double tol = ABS_TOL) {
    for (size_t i = 0; i < N; ++i) {
        if (!approx_equal(a[i], b[i], tol)) {
            return false;
        }
    }

    return true;
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
