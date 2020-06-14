#pragma once

#include "IntegralFrechet/IntegralFrechet.h"
#include "Frechet/frechet_light.h"
#include "DTW/dtw.h"

distance_t dtw(Curve curve_1, Curve curve_2) {
    return DTW(curve_1, curve_2).cost();
}

distance_t integral_frechet(Curve curve_1, Curve curve_2) {
    // auto start = std::chrono::high_resolution_clock::now();
    auto dist = IntegralFrechet(
        curve_1, curve_2, ParamMetric::L1, 5, nullptr
    ).compute_matching()
    .cost;
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
    return dist;
}

distance_t integral_frechet_fast(Curve curve_1, Curve curve_2) {
    const auto result_alt = IntegralFrechet(curve_1.simplify(true), curve_2, ParamMetric::L1, 1, nullptr).compute_matching();
    const auto band = MatchingBand(curve_1, curve_2, result_alt.matching, 1);
    const auto result = IntegralFrechet(curve_1, curve_1, ParamMetric::L1, 1, &band).compute_matching();
    return result.cost;
}

distance_t frechet(Curve curve_1, Curve curve_2) {
    FrechetLight fl;
    return fl.calcDistance(curve_1, curve_2);
}

distance_t average_frechet(Curve curve_1, Curve curve_2) {
    return (1 / (curve_1.curve_length() + curve_2.curve_length())) * integral_frechet(curve_1, curve_2);
};

distance_t average_frechet_fast(Curve curve_1, Curve curve_2) {
    return (1 / (curve_1.curve_length() + curve_2.curve_length())) * integral_frechet_fast(curve_1, curve_2);
}

