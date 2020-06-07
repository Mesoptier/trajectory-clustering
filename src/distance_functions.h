#ifndef DISTANCE_FUNCTIONS
#define DISTANCE_FUNCTIONS

#include "Curve.h"

namespace df {
    distance_t dtw(const Curve& curve_1, const Curve& curve_2);

    distance_t integral_frechet(const Curve& curve_1, const Curve& curve_2);

    distance_t integral_frechet_fast(const Curve& curve_1,
        const Curve& curve_2);

    distance_t frechet(const Curve& curve_1, const Curve& curve_2);

    distance_t average_frechet(const Curve& curve_1, const Curve& curve_2);

    distance_t average_frechet_fast(const Curve& curve_1, const Curve& curve_2);

    bool dtw_lt(const Curve& curve_1, const Curve& curve_2, distance_t delta);

    bool integral_frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta);

    bool integral_frechet_fast_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta);

    bool frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta);

    bool average_frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta);

    bool average_frechet_fast_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta);
}
#endif
