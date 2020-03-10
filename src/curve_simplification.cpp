#include "curve_simplification.h"
#include "defs.h"
#include <limits>
#include "IntegralFrechet/IntegralFrechet.h"

namespace
{

distance_t compute_integral_frechet_distance(Curve curve_1, Curve curve_2) {
    std::cout << "computing distance...simplification \n";
    return IntegralFrechet(
        curve_1, curve_2, ParamMetric::LInfinity_NoShortcuts, 1, nullptr
    ).compute_matching()
    .cost;
}

Curve simplify(Curve const& curve, distance_t distance)
{
	Curve simplified_curve({curve.front()});
	Curve prefix_curve({curve.front()});

	// FrechetLight frechet_light;
	for (PointID id = 1; id < curve.size()-1; ++id) {
		auto const& point = curve[id];
		prefix_curve.push_back(point);
		auto line_segment = Curve({prefix_curve.front(), prefix_curve.back()});
		bool less_than = distance < compute_integral_frechet_distance(line_segment, prefix_curve);
		// bool less_than = false;
		// bool less_than = frechet_light.lessThanWithFilters(distance, line_segment, prefix_curve);
		if (!less_than) {
			simplified_curve.push_back(point);
			prefix_curve = Curve({point});
		}
	}
	simplified_curve.push_back(curve.back());

	return simplified_curve;
}

} // end anonymous namespace

// triest to find an l-simplification with a small distance
Curve simplify(Curve const& curve, int l)
{
	static constexpr distance_t epsilon = 1e-8;

	assert(l >= 2);
	if ((int)curve.size() < l) { return curve; }

	distance_t min = 0.;
	// distance_t max = curve.getUpperBoundDistance(curve);
	distance_t max = MAXFLOAT;

	Curve simplified_curve = curve;
	while (max-min > epsilon) {
		auto split = (max + min)/2.;
		simplified_curve = simplify(curve, split);
		if ((int)simplified_curve.size() <= l) {
			max = split;
		}
		else {
			min = split;
		}
	}

	return simplify(curve, max);
}
