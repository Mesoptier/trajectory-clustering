#include "distance_functions.h"

#include "DTW/dtw.h"
#include "Frechet/frechet_light.h"
#include "IntegralFrechet/IntegralFrechet.h"

distance_t df::dtw(const Curve& curve_1, const Curve& curve_2) {
    return DTW(curve_1, curve_2).cost();
}

distance_t df::integral_frechet(const Curve& curve_1, const Curve& curve_2) {
    return IntegralFrechet(curve_1, curve_2, ParamMetric::L1, 1, nullptr)
    	.compute_matching().cost;
}

distance_t df::integral_frechet_fast(const Curve& curve_1,
		const Curve& curve_2) {
    const auto result_alt = IntegralFrechet(curve_1.simplify().curve, curve_2,
    	ParamMetric::L1, 1, nullptr).compute_matching();
    const auto band = MatchingBand(curve_1, curve_2, result_alt.matching, 1);
    return IntegralFrechet(curve_1, curve_1, ParamMetric::L1, 1,&band)
    	.compute_matching().cost;
}

distance_t df::frechet(const Curve& curve_1, const Curve& curve_2) {
    FrechetLight fl;
    return fl.calcDistance(curve_1, curve_2);
}

distance_t df::average_frechet(const Curve& curve_1, const Curve& curve_2) {
    return integral_frechet(curve_1, curve_2) /
    	(curve_1.curve_length() + curve_2.curve_length());
}

distance_t df::average_frechet_fast(const Curve& curve_1,
		const Curve& curve_2) {
    return integral_frechet_fast(curve_1, curve_2) /
    	(curve_1.curve_length() + curve_2.curve_length());
}

bool df::dtw_lt(const Curve& curve_1, const Curve& curve_2, distance_t delta) {
	return dtw(curve_1, curve_2) <= delta;
}

bool df::integral_frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta) {
	return integral_frechet(curve_1, curve_2) <= delta;
}

bool df::integral_frechet_fast_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta) {
	return integral_frechet_fast(curve_1, curve_2) <= delta;
}

bool df::frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta) {
	FrechetLight fl;
	return fl.lessThanWithFilters(delta, curve_1, curve_2);
}

bool df::average_frechet_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta) {
	return average_frechet(curve_1, curve_2) <= delta;
}

bool df::average_frechet_fast_lt(const Curve& curve_1, const Curve& curve_2,
        distance_t delta) {
	return average_frechet_fast(curve_1, curve_2) <= delta;
}
