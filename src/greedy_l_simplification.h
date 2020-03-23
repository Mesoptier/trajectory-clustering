#pragma once

#include "greedy_simplification.h"

Curve greedy_l_simplification(Curve const& curve, std::size_t l) {

    distance_t epsilon = 1e-8;
    distance_t min = 0;
    distance_t max = IntegralFrechet(curve, Curve("", {curve.get_points()[0], curve.get_points()[curve.get_points().size()-1]}), ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
    .compute_matching().cost;

    std::cout << "max: " << max << "\n";

    Curve simplified_curve = curve;
	while (max-min > epsilon) {
        std::cout << max-min << "\n";
		auto split = (max + min)/2.;
		simplified_curve = greedy_simplification(curve, split);
		if (simplified_curve.size() <= l) {
			max = split;
		}
		else {
			min = split;
		}
    }
    return simplified_curve;
}
