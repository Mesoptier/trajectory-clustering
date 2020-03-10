#pragma once

#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"

distance_t compute_integral_frechet_distance(Curve curve_1, Curve curve_2) {
	/*const auto result_alt = IntegralFrechet(curve_1.simplify(true), curve_2, ParamMetric::LInfinity_NoShortcuts, 10, nullptr).compute_matching();
    const auto band = MatchingBand(curve_1, curve_2, result_alt.matching, 1);
	const auto result = IntegralFrechet(curve_1, curve_1, ParamMetric::LInfinity_NoShortcuts, 1, &band).compute_matching();


	std::cout << result.cost << "\n";
	return result.cost;*/
	
	distance_t cost = IntegralFrechet(
        curve_1, curve_2, ParamMetric::LInfinity_NoShortcuts, 1, nullptr
    ).compute_matching()
    .cost;

	std::cout << cost << "\n";

    return cost;
}
