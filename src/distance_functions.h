#pragma once

#include "IntegralFrechet/IntegralFrechet.h"


distance_t integral_frechet(Curve curve_1, Curve curve_2) {

    return IntegralFrechet(
        curve_1, curve_2, ParamMetric::LInfinity_NoShortcuts, 1, nullptr
    ).compute_matching()
    .cost;

}

distance_t frechet(Curve curve_1, Curve curve_2) {
    return 0;
}

distance_t dtw(Curve curve_1, Curve curve_2) {
    return 0;
}

