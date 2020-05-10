#include "curve_simplification.h"
#include "defs.h"
#include <limits>
#include "IntegralFrechet/IntegralFrechet.h"


namespace
{

Curve simplify(Curve const& curve, distance_t distance, distance_t(*dist_func)(Curve, Curve))
{
	std::cout << "using curve simplification...\n";
	Curve simplified_curve({curve.front()});
	Curve prefix_curve({curve.front()});

	FrechetLight fl;

	for (PointID id = 1; id < curve.size()-1; ++id) {
		auto const& point = curve[id];
		prefix_curve.push_back(point);
		auto line_segment = Curve({prefix_curve.front(), prefix_curve.back()});
		bool less_than = dist_func(line_segment, prefix_curve) < distance;
		// bool less_than = fl.lessThanWithFilters(distance, line_segment, prefix_curve);
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
Curve simplify(Curve const& curve, std::size_t l, distance_t(*dist_func)(Curve, Curve))
{
	static constexpr distance_t epsilon = 1e-8;
	assert(l >= 2);
	if (curve.size() < l) { return curve; }

	distance_t min = 0.;
	distance_t max = curve.getUpperBoundDistance(curve);

	Curve simplified_curve = curve;
	while (max-min > epsilon) {
		auto split = (max + min)/2.;
		simplified_curve = simplify(curve, split, dist_func);
		if (simplified_curve.size() <= l) {
			max = split;
		}
		else {
			min = split;
		}
	}

	return simplify(curve, max, dist_func);
}
