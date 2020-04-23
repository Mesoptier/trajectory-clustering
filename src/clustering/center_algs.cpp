#include "center_algs.h"

//#include "curve_simplification.h"
#include "../Frechet/frechet_light.h"
#include "../Frechet/frechet_matching.h"
#include "../geom.h"
#include "../IntegralFrechet/IntegralFrechet.h"

#include <limits>

namespace
{

Points get_uniformly_spaced_points(Points param_space_path, Curve curve_1, Curve curve_2) {

	// std::cout << "uniform points...\n";
	// std::cout << param_space_path.size() << "\n";
	// std::cout << "curve_1 size... " << curve_1.get_points().size() << "\n"; 
	//sequence of points on curve_2
	Points points = Points();

	distance_t curve_1_length = curve_1.curve_length();
	distance_t interval_size = curve_1_length / (curve_1.get_points().size() - 1);

	std::vector<distance_t> x_coords = std::vector<distance_t>();
	for (auto p: param_space_path) {
		x_coords.push_back(p.x);
	}
	points.push_back(curve_2[0]);
	for (int i = 1; i < curve_1.get_points().size(); ++i) {
		auto it =  std::lower_bound(x_coords.begin(), x_coords.end(), interval_size * i);
		auto prev_index = static_cast<std::size_t>(std::distance(x_coords.begin(), it)) - 1;
		auto next_index = prev_index+1;
		
		if (prev_index >= param_space_path.size()) {
			std::cout << "yes\n";
			prev_index = param_space_path.size() - 2;
			next_index = prev_index + 1;
		}

		Point p = param_space_path[prev_index];
		Point q = param_space_path[next_index];

		// std::cout << "p: " << p << "\n";
		// std::cout << "q: " << q << "\n";
		// std::cout << interval_size * (i+1) << "\n";
		// std::cout << "prev: " << prev_index << "\n";
		// std::cout << "next: " << next_index << "\n";
		Line line = Line::fromTwoPoints(p, q);
		distance_t distance_along_c2;
		if (line.isVertical()) {
			distance_along_c2 = p.y;
		} else {
			distance_along_c2 = line.getY(interval_size * i);
			assert(distance_along_c2 > 0);
		}

		if (distance_along_c2 > curve_2.curve_length()) {
			distance_along_c2 = curve_2.curve_length();
		}
		CPoint cpoint = curve_2.get_cpoint_after(distance_along_c2);
		Point new_point = curve_2.interpolate_at(cpoint);

		// if (points.size() > 0)
		// 	assert(!approx_equal(points.back(), new_point));
		
		points.push_back(
			new_point
		);
	}

	return points;

}

Points matching_to_points(Points param_space_path, Curve curve_1, Curve curve_2) {
	std::vector<distance_t> curve_2_lengths = std::vector<distance_t>();
	curve_2_lengths.push_back(0);

	std::vector<distance_t> x_coords = std::vector<distance_t>();
	for (auto p: param_space_path) {
		x_coords.push_back(p.x);
	}

	for (std::size_t i = 1; i < curve_1.size(); ++i) {
		distance_t arc_length = curve_1.curve_length(i);

		auto it =  std::lower_bound(x_coords.begin(), x_coords.end(), arc_length);
		auto index = static_cast<std::size_t>(std::distance(x_coords.begin(), it));

		curve_2_lengths.push_back(param_space_path[index].y);
	}


	// Find the points on curve_2 corresponding to the sequence of arc lengths
	Points points = Points();

	for (auto length: curve_2_lengths) {
		if (length > curve_2.curve_length()) {
			// std::cout << length - curve_2.curve_length() << "\n";
			length = curve_2.curve_length();
		}
		assert(length <= curve_2.curve_length());
		CPoint cpoint = curve_2.get_cpoint_after(length);
		Point point = curve_2.interpolate_at(cpoint);
		points.push_back(point);
	}

	return points;
}


Point mean_of_points(Points points) {
	distance_t x_mean = 0;
	distance_t y_mean = 0;

	for (auto p: points) {
		x_mean += p.x;
		y_mean += p.y;
	}

	x_mean /= points.size();
	y_mean /= points.size();

	return Point(x_mean, y_mean);
}

bool calcKXCenters(Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve))
{
	bool found_new_center = false;

	// compute cluster costs in case they haven't been computed yet
	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

	for (auto& cluster: clustering) {
		for (CurveID curve_id1: cluster.curve_ids) {
			// auto simplified_curve = simplify(curves[curve_id1], l);
			auto simplified_curve = curves[curve_id1].naive_l_simplification(l);
			auto dist = calcC2CDist(curves, simplified_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (dist < cluster.cost) {
				cluster.center_curve = simplified_curve;
				cluster.cost = dist;
				found_new_center = true;
			}
		}
	}

	return found_new_center;
}

} // end anonymous namespace

std::string toString(CenterAlg center_alg) {
	switch(center_alg) {
	case CenterAlg::kMedian: return "kMedian";
	case CenterAlg::kMeans: return "kMeans";
	case CenterAlg::kCenter: return "kCenter";
	case CenterAlg::fCenter: return "fCenter";
	case CenterAlg::fMean: return "fMean";
	}
	ERROR("Unknown center_alg.");
}

distance_t calcC2CDist(
	Curves const& curves, Curve const& center_curve, CurveIDs const& curve_ids, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve))
{
	// FrechetLight frechet_light;

	distance_t dist = 0;
	for (auto curve_id: curve_ids) {
		// auto curve_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
		// auto curve_dist = IntegralFrechet(center_curve, curves[curve_id], ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
		// .compute_matching()
		// .cost;
		auto curve_dist = dist_func(center_curve, curves[curve_id]);

		switch (c2c_dist) {
		case C2CDist::Median:
			dist += curve_dist;
			break;
		case C2CDist::Mean:
			dist += curve_dist*curve_dist;
			break;
		case C2CDist::Max:
			dist = std::max(dist, curve_dist);
			break;
		}
	}

	return dist;
}

bool computerCenters(Curves const& curves, Clustering& clustering, int l, CenterAlg center_alg, distance_t(*dist_func)(Curve, Curve))
{
	switch (center_alg) {
	case CenterAlg::kMedian:
		return calcKMedianCenters(curves, clustering, l, dist_func);
	case CenterAlg::kMeans:
		return calcKMeansCenters(curves, clustering, l, dist_func);
	case CenterAlg::kCenter:
		return calcKCenterCenters(curves, clustering, l, dist_func);
	case CenterAlg::fCenter:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Max, CenterCurveUpdateMethod::frechetCentering);
	case CenterAlg::fMean:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::frechetMean);
	case CenterAlg::dtwMean:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::dtwMean);
	case CenterAlg::avFCenter:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::avFCenter);
	case CenterAlg::newCenterUpdate:
		return newCenterUpdate(curves, clustering, dist_func, C2CDist::Median);
	}

	ERROR("No matching center_alg enum passed.");
}

bool calcKMedianCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve))
{
	return calcKXCenters(curves, clustering, l, C2CDist::Median, dist_func);
}

bool calcKMeansCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve))
{
	return calcKXCenters(curves, clustering, l, C2CDist::Mean, dist_func);
}

bool calcKCenterCenters(Curves const& curves, Clustering& clustering, int l, distance_t(*dist_func)(Curve, Curve))
{
	return calcKXCenters(curves, clustering, l, C2CDist::Max, dist_func);
}

bool frechetCentering (Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve)) {
	return calcFSACenters(curves, clustering, l, dist_func, c2c_dist, CenterCurveUpdateMethod::frechetCentering);
}

bool frerchetMean(Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist, distance_t(*dist_func)(Curve, Curve)) {
	return calcFSACenters(curves, clustering, l, dist_func, c2c_dist, CenterCurveUpdateMethod::frechetMean);
}

// TODO: There is some unnecessary pushing around of data here. Fix that to increase performance.
bool calcFSACenters(Curves const& curves, Clustering& clustering, int /* l */, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist, CenterCurveUpdateMethod method)
{
	bool found_new_center = false;
	// FrechetLight frechet_light;

	// compute cluster costs in case they haven't been computed yet
	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

	for (auto& cluster: clustering) {
		std::vector<Points> matchings;
		auto const& center_curve = cluster.center_curve;
		Curve new_center_curve;

		for (auto curve_id: cluster.curve_ids) {
			auto const& curve = curves[curve_id];
			// auto matching = calcMatching(cluster.center_curve, curve);
			// auto matching = IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
			// .compute_matching()
			// .matching;

			// auto matching = method == frechetCentering ? 
			// calcMatching(cluster.center_curve, curve) :
			// IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
			// .compute_matching()
			// .matching;

			// Points matching_points = matching_to_points(matching, cluster.center_curve, curve);
			// matchings.push_back(std::move(matching_points));

			switch (method) {
				case CenterCurveUpdateMethod::frechetCentering:
					matchings.push_back(calcMatching(cluster.center_curve, curve));
					break;
				case CenterCurveUpdateMethod::avFCenter:
					matchings.push_back(
						matching_to_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
							.compute_matching()
							.matching,
							cluster.center_curve,
							curve
						)
					);
					break;
				case CenterCurveUpdateMethod::frechetMean:
					matchings.push_back(
						matching_to_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
							.compute_matching()
							.matching,
							cluster.center_curve,
							curve
						)
					);
					break;
				case CenterCurveUpdateMethod::dtwMean:
					auto dtw_matching = DTW(cluster.center_curve, curve).matching();
					Points matching = Points();
					int index = 0;
					for (int i = 0; i < cluster.center_curve.get_points().size(); ++i) {
						while (i != dtw_matching[index].first) {
							index++;
						}
						matching.push_back(curve[index]);
					}
					break;
			}
		}

		for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
			Points matching_points;
			for (auto const& matching: matchings) {
				matching_points.push_back(matching[point_id]);
			}

			Point new_point;

			switch (method) {
				case CenterCurveUpdateMethod::frechetCentering:
					new_point = calcMinEnclosingCircle(matching_points).center;
					break;
				case CenterCurveUpdateMethod::avFCenter:
					new_point = calcMinEnclosingCircle(matching_points).center;
					break;
				case CenterCurveUpdateMethod::frechetMean:
					new_point = mean_of_points(matching_points);
					break;
				case CenterCurveUpdateMethod::dtwMean:
					new_point = mean_of_points(matching_points);
					break;
			}

			if (new_center_curve.get_points().empty() || !approx_equal(new_point, new_center_curve.get_points().back())) {
				new_center_curve.push_back(new_point);
			}
		}

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (new_dist < cluster.cost) {
				cluster.center_curve = std::move(new_center_curve);
				cluster.cost = new_dist;
				found_new_center = true;
			}
		} 
	}
	
	if (found_new_center) {
		std::cout << "found new center\n";
	}
	else {
		std::cout << "no new center... :( \n";
	}
	return found_new_center;
}

bool newCenterUpdate(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
	bool found_new_center = false;

	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

		for (auto& cluster: clustering) {
		std::vector<Points> matchings;
		auto const& center_curve = cluster.center_curve;
		Curve new_center_curve;

		for (auto curve_id: cluster.curve_ids) {
			auto const& curve = curves[curve_id];

			matchings.push_back(
				get_uniformly_spaced_points(
					IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
					.compute_matching()
					.matching,
					cluster.center_curve,
					curve
				)
			);

		}

		for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
			Points matching_points;
			for (auto const& matching: matchings) {
				matching_points.push_back(matching[point_id]);
			}

			Point new_point;
			new_point = mean_of_points(matching_points);


			if (new_center_curve.get_points().empty() || !approx_equal(new_point, new_center_curve.get_points().back())) {
				new_center_curve.push_back(new_point);
			}
		}

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (new_dist < cluster.cost) {
				cluster.center_curve = std::move(new_center_curve);
				cluster.cost = new_dist;
				found_new_center = true;
			}
		} 
	}

		
	if (found_new_center) {
		std::cout << "found new center\n";
	}
	else {
		std::cout << "no new center... :( \n";
	}
	return found_new_center;


}

Points matching_of_vertices(Curve curve_1, Curve curve_2) {
		auto matching = IntegralFrechet(curve_1, curve_2, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
		.compute_matching().matching;
		return matching_to_points(matching, curve_1, curve_2);
	}

