#include "center_algs.h"
#include "wedge_method.h"
#include "regression_method.h"
//#include "curve_simplification.h"
#include "../Frechet/frechet_light.h"
#include "../Frechet/frechet_matching.h"
#include "../geom.h"
#include "../IntegralFrechet/IntegralFrechet.h"
#include <limits>

namespace
{

std::vector<distance_t> get_redistributed_distances(Curve curve, 
		double lower_threshold=1./10, double upper_threshold=1./5) {

	Points new_points = Points();

	//indices of new points
	std::vector<size_t> indices = std::vector<size_t>();
	
	int point_budget = 0;

	// First pass over curve
	// Remove redundant points
	for (int i = 0; i < curve.get_points().size(); ++i) {
		if (i > 0 && i < curve.get_points().size() - 1) {

			Point a = new_points.back();
			Point b = curve[i];
			Point c = curve[i+1];

			Point dir_1 = Point(b.x - a.x, b.y - a.y);
			Point dir_2 = Point(c.x - b.x, c.y - b.y);

			distance_t theta = acute_angle(dir_1, dir_2);
			// std::cout << theta << " >= " << lower_threshold << "\n";
			if (theta >= lower_threshold) {
				if (!approx_equal(new_points.back(), curve[i])) {
					new_points.push_back(curve[i]);
					indices.push_back(i);
					continue;
				}
			} else {
				// std::cout << "skipped a point\n";
			}

			point_budget++;

		} else {
			if (i == 0 || !approx_equal(new_points.back(), curve[i])) {
				new_points.push_back(curve[i]);
				indices.push_back(i);
			}
		}
	}

	// Second pass over the curve, add extra
	// points to sections with high curvature
	Points final_points = Points();
	std::vector<distance_t> final_distances = std::vector<distance_t>();
	size_t last_index = 0;

	int i = 0;

	while (i < new_points.size()) {
		if (i > 0 && i < new_points.size() - 1) {
			
			Point a = final_points.back();
			Point b = new_points[i];
			Point c = new_points[i+1];

			Point dir_1 = Point(b.x - a.x, b.y - a.y);
			Point dir_2 = Point(c.x - b.x, c.y - c.y);

			distance_t theta = acute_angle(dir_1, dir_2);

			if (theta > upper_threshold && point_budget > 0) {
				std::cout << "new point...\n";
				distance_t dist_1 = a.dist(b);
				distance_t dist_2 = b.dist(c);

				assert(dist_1 > 0);
				assert(dist_2 > 0);

				if (dist_1 > dist_2) {
					final_distances.push_back(final_distances.back() + 9* dist_1 / 10);
					final_distances.push_back(final_distances.back() + dist_1 / 10);
					final_distances.push_back(final_distances.back() + dist_2);
				} else {
					final_distances.push_back(final_distances.back() + dist_1);
					final_distances.push_back(final_distances.back() + dist_2 / 10);
					final_distances.push_back(final_distances.back() + 9 * dist_2 / 10);
				}

				final_points.push_back(c);
				point_budget--;

				// Curve temp_curve = Curve({a, b, c});
				// distance_t total_length = temp_curve.curve_length();
				// distance_t third = total_length / 3;

				// CPoint cpoint_1 = temp_curve.get_cpoint_after(third);
				// CPoint cpoint_2 = temp_curve.get_cpoint_after(2 * third);

				// Point new_point_1 = temp_curve.interpolate_at(cpoint_1);
				// Point new_point_2 = temp_curve.interpolate_at(cpoint_2);

				// final_points.push_back(new_point_1);
				// final_points.push_back(new_point_2);
				// final_points.push_back(c);

				// distance_t dist_a = curve.curve_length(last_index);
				// distance_t new_dist_1 = curve.curve_length(
				// 	curve.get_cpoint_after(
				// 		curve.curve_length(last_index) + third
				// 	)
				// );
				// distance_t new_dist_2 = curve.curve_length(
				// 	curve.get_cpoint_after(
				// 		curve.curve_length(last_index) + 2 * third
				// 	)
				// );
				// distance_t dist_c = curve.curve_length(indices[i+1]);


				// std::cout << "extra points!!";
				// final_distances.push_back(new_dist_1);
				// final_distances.push_back(new_dist_2);
				// final_distances.push_back(dist_c);


				// last_index = indices[i+1];
				// point_budget--;
			} else {
				final_points.push_back(b);
				final_points.push_back(c);
				last_index = indices[i+1];

				final_distances.push_back(curve.curve_length(indices[i]));
				final_distances.push_back(curve.curve_length(indices[i+1]));
			}

			i += 2;

		} else {
			if (final_distances.size() > 0) {
				assert(!approx_equal(final_distances.back(), curve.curve_length(indices[i])));
			}
			final_points.push_back(new_points[i]);
			final_distances.push_back(curve.curve_length(indices[i]));
			i += 1;
		}
	}

	// for (auto dist: final_distances) {
	// 	std::cout << dist << "\n";
	// }
	// std::cout << Points(curve.get_points()) << "\n";
	// std::cout << "length: " << curve.curve_length();
	return final_distances;
}

// Finds points on curve_2 alined with
// the points on another curve defined by distances.
// The param_space_path argument is a sequence of 
// points defining the matching between the two curves.

Points get_images_of_points(
	Points param_space_path, std::vector<distance_t> distances,  Curve curve_2
) {
	Points points = Points();

	std::vector<distance_t> x_coords = std::vector<distance_t>();
	for (auto p: param_space_path) {
		x_coords.push_back(p.x);
	}
	points.push_back(curve_2[0]);


	int last_index = 0;

	for (int i = 0; i < distances.size(); ++i) {

		distance_t distance = distances[i];
		
		if (distance == 0)
			continue;

		// while (last_index <= x_coords.size() - 2 && x_coords[last_index + 1] <= distance) {
		// 	++last_index;
		// }

		auto it = std::lower_bound(x_coords.begin(), x_coords.end(), distance);
		last_index = static_cast<std::size_t>(std::distance(x_coords.begin(), it));

		int prev = last_index;
		int next = last_index + 1;

		if (prev >= param_space_path.size()) {
			prev = param_space_path.size() - 2;
			next = prev + 1;
		}

		Point p = param_space_path[prev];
		Point q = param_space_path[next];

		Line line = Line::fromTwoPoints(p, q);
		distance_t distance_along_c2;
		if (line.isVertical()) {
			distance_along_c2 = p.y;
		} else {
			distance_along_c2 = line.getY(distances[i]);
			// std::cout << distance_along_c2 << ", " << curve_2.curve_length() << "\n";
			assert(distance_along_c2 >= 0);
			// assert(distance_along_c2 <= param_space_path.back().y);
		}

		if (distance_along_c2 > curve_2.curve_length()) {
			distance_along_c2 = curve_2.curve_length();
		}
		CPoint cpoint = curve_2.get_cpoint_after(distance_along_c2);
		Point new_point = curve_2.interpolate_at(cpoint);

		// if (i == distances.size() - 1) {
		// 	if (!approx_equal(new_point, curve_2.get_points().back())) {
		// 		std::cout << "uh oh...\n";
		// 	}
		// 	// assert(approx_equal(new_point, curve_2.get_points().back()));
		// }

		points.push_back(
			new_point
		);
	}

	return points;
}

std::vector<distance_t> get_uniform_distances_along_curve(Curve curve) {

	distance_t distance = curve.curve_length() / curve.size();

	return std::vector<distance_t>(curve.size() - 1, distance);
}

Points get_uniformly_spaced_points(Points param_space_path, Curve curve_1, Curve curve_2) {

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
			prev_index = param_space_path.size() - 2;
			next_index = prev_index + 1;
		}

		Point p = param_space_path[prev_index];
		Point q = param_space_path[next_index];

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
		
		points.push_back(
			new_point
		);
	}

	return points;
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
	case CenterAlg::newCenterUpdate: return "newCenterUpdate";
	case CenterAlg::newCenterUpdate2: return "newCenterUpdate2";
	case CenterAlg::ensembleMethod1: return "ensembleMethod1";
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
	case CenterAlg::dba:
		return dba(curves, clustering, dist_func, C2CDist::Mean);
	case CenterAlg::avFCenter:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::avFCenter);
	case CenterAlg::newCenterUpdate:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::newCenterUpdate);
	case CenterAlg::newCenterUpdate2:
		return calcFSACenters(curves, clustering, l, dist_func, C2CDist::Median, CenterCurveUpdateMethod::newCenterUpdate2);
	case CenterAlg::naiveCenterUpdate:
		return naiveCenterUpdate(curves, clustering, dist_func, C2CDist::Median);
	case CenterAlg::ensembleMethod1:
		return ensembleMethod1(curves, clustering, dist_func, C2CDist::Median);
	case CenterAlg::cdba:
		return cdba(curves, clustering, dist_func, C2CDist::Median);
	case CenterAlg::wedge:
		return wedge_method(curves, clustering, dist_func, C2CDist::Median);
	case CenterAlg::regression:
		return regression_method(curves, clustering, dist_func, C2CDist::Median);
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
						get_images_of_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::L1, 1, nullptr)
							.compute_matching()
							.matching,
							cluster.center_curve.get_prefix_length_vector(),
							curve
						)
					);
					break;
				case CenterCurveUpdateMethod::frechetMean:
					matchings.push_back(
						get_images_of_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::L1, 1, nullptr)
							.compute_matching()
							.matching,
							cluster.center_curve.get_prefix_length_vector(),
							curve
						)
					);
					break;
				case CenterCurveUpdateMethod::dtwMean: {
					auto dtw_matching = DTW(cluster.center_curve, curve).matching();
					Points matching = Points();
					int index = 0;
					for (int i = 0; i < cluster.center_curve.get_points().size(); ++i) {
						while (i != dtw_matching[index].first) {
							++index;
						}
						while (i == dtw_matching[index].first && index < dtw_matching.size()) {
							matching.push_back(curve[dtw_matching[index].second]);
							++index;
						}
					}
					matchings.push_back(matching);
					break;
				}
				case CenterCurveUpdateMethod::newCenterUpdate:
					matchings.push_back(
						get_images_of_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
							.compute_matching()
							.matching,
							get_uniform_distances_along_curve(cluster.center_curve),
							curve
						)
					);
					break;
				case CenterCurveUpdateMethod::newCenterUpdate2:
					matchings.push_back(
						get_images_of_points(
							IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
							.compute_matching()
							.matching,
							get_redistributed_distances(cluster.center_curve),
							curve
						)
					);
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
				case CenterCurveUpdateMethod::newCenterUpdate:
					new_point = mean_of_points(matching_points);
					break;
				case CenterCurveUpdateMethod::newCenterUpdate2:
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
				std::cout << "new cluster center\n";
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

bool dba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
	bool found_new_center = false;

	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

	for (auto& cluster: clustering) {
		const auto& center_curve = cluster.center_curve;
		Curve new_center_curve;
		std::vector<std::vector<Points>> matchings = std::vector<std::vector<Points>>();

		for (auto& curve_id: cluster.curve_ids) {
			Curve curve = curves[curve_id];
			auto dtw_matching = DTW(cluster.center_curve, curve).matching();
			std::vector<Points> matching = std::vector<Points>();

			int matching_index = 0;
			for (int i = 0; i < cluster.center_curve.size(); ++i) {
				matching.push_back({});
				while (matching_index < dtw_matching.size() && dtw_matching[matching_index].first == i) {
					matching.back().push_back(curve[dtw_matching[matching_index].second]);
					++matching_index;
				}
			}

			matchings.push_back(matching);
		}

		for (int i = 0; i < cluster.center_curve.size(); ++i) {
			Points points_to_average = Points();

			for (auto& matching: matchings) {
				// points_to_average.push_back(matching[i][0]);
				for (auto& point: matching[i]) {
					points_to_average.push_back(point);
				}
			}

			Point new_point = mean_of_points(points_to_average);

			if (new_center_curve.size() == 0 || !approx_equal(new_point, new_center_curve.back())) {
				new_center_curve.push_back(
					new_point
				);
			}

		}

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (new_dist < cluster.cost) {
				std::cout << "new cluster center\n";
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

std::pair<distance_t, distance_t> get_y_range(Points& param_space_path, distance_t distance) {

	std::vector<distance_t> x_coords = std::vector<distance_t>();

	for (auto& point: param_space_path) {
		x_coords.push_back(point.x);
	}

	if (distance > param_space_path.back().x)
		distance = param_space_path.back().x;

	int index = 0;

	auto it = std::lower_bound(x_coords.begin(), x_coords.end(), distance);
	index = static_cast<std::size_t>(std::distance(x_coords.begin(), it));

	if (x_coords[index] > distance)
		index--;

	// if (index < 0)
	// 	++index;

	// while (index <= x_coords.size() - 2 && (x_coords[index + 1] <= distance)) {
	// 	index++;
	// }

	if (index == x_coords.size() - 1) {
		return {param_space_path.back().y, param_space_path.back().y};
	}

	Point p = param_space_path[index];
	Point q = param_space_path[index + 1];

	if (approx_equal(p, q)) {
		std::cout << "...\n";
	}

	Line line = Line::fromTwoPoints(p, q);

	if (!line.isVertical()) {
		distance_t num = line.getY(distance);

		if (isnan(num)) {
			std::cout << "num is nan...\n";
		}

		if (approx_equal(num, 0.))
			num = 0.;

		if (num  < 0)
			std::cout << "this is weird...\n";

		return {num, num};
	}

	int next_index = index + 1;

	while (next_index < param_space_path.size() && param_space_path[next_index].x == param_space_path[index].x)
		++next_index;
	

	if (isnan(param_space_path[index].y))
		std::cout << "uh-oh spaghetti-ohs\n";
	
	return {param_space_path[index].y, param_space_path[next_index - 1].y};
}

bool cdba(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
	bool found_new_center = false;

	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

	for (auto& cluster: clustering) {
		const auto& center_curve = cluster.center_curve;
		Curve new_center_curve;
		std::vector<std::vector<Points>> matchings = std::vector<std::vector<Points>>();

		for (auto& curve_id: cluster.curve_ids) {
			Curve curve = curves[curve_id];
			Points param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 1, nullptr)
			.compute_matching()
			.matching;

			std::vector<Points> matching = std::vector<Points>();

			for (int i = 0; i < center_curve.size(); ++i) {
				matching.push_back({});
				distance_t dist = center_curve.curve_length(i);
				std::pair<distance_t, distance_t> y_range = get_y_range(param_space_path, dist);

				// if (i == 0) {
				// 	assert(approx_equal(curve.interpolate_at(curve.get_cpoint_after(y_range.first)), {0, 0}));
				// }

				if (approx_equal(y_range.first, y_range.second)) {
					matching.back().push_back(curve.interpolate_at(curve.get_cpoint_after(y_range.first)));
				}
				else {
					distance_t length = y_range.second - y_range.first;
					distance_t ratio_to_curve_length = length / curve.curve_length();
					int number_of_samples = 2 * curve.size() * ratio_to_curve_length;

					if (number_of_samples == 0)
						number_of_samples = 1;

					for (int j = 0; j < number_of_samples; ++j) {
						distance_t dist_along_curve = y_range.first + j * length / number_of_samples;
						matching.back().push_back(curve.interpolate_at(curve.get_cpoint_after(dist_along_curve)));
					}
				}

				if (matching.back().size() == 0)
					std::cout << "woah...\n";
			}

			matchings.push_back(matching);
		}

		for (int i = 0; i < cluster.center_curve.size(); ++i) {
			Points points_to_average = Points();

			for (auto& matching: matchings) {
				// points_to_average.push_back(matching[i][0]);
				for (auto& point: matching[i]) {
					points_to_average.push_back(point);
				}
			}

			Point new_point = mean_of_points(points_to_average);

			if (new_center_curve.size() == 0 || !approx_equal(new_point, new_center_curve.back())) {
				new_center_curve.push_back(
					new_point
				);
			}

		}

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (new_dist < cluster.cost) {
				// std::cout << "new cluster center\n";
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

bool naiveCenterUpdate(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
	bool found_new_center = false;

	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

	for (auto& cluster: clustering) {
		auto& center_curve = cluster.center_curve;
		int l = center_curve.get_points().size();
		std::vector<Points> matchings = std::vector<Points>();
		for (int i = 0; i < l; ++i) {
			matchings.push_back(Points());
		}
	}

	for (auto& cluster: clustering) {
		auto& center_curve = cluster.center_curve;
		int l = center_curve.get_points().size();
		std::vector<Points> matchings = std::vector<Points>();
		for (int i = 0; i < l; ++i) {
			matchings.push_back(Points());
		}


		for (auto& curve_id: cluster.curve_ids) {
			Curve curve = curves[curve_id];
			distance_t dist = curve.curve_length() / (l - 1);

			for (int i = 0; i < l; ++i) {
				matchings[i].push_back(curve.interpolate_at(curve.get_cpoint_after(dist * i)));
			}
		}
		std::cout << "finished first loop\n";
		Points new_points = Points();

		for (size_t i = 0; i < matchings.size(); ++i) {
			if (matchings[i].size() > 0) {
				Point point = mean_of_points(matchings[i]);
				// if (new_points.size() > 0 && !approx_equal(new_points.back(), point))
				new_points.push_back(point);
			}
		}

		Curve new_center_curve = Curve(new_points);

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
			if (new_dist < cluster.cost) {
				cluster.center_curve = std::move(new_center_curve);
				cluster.cost = new_dist;
				found_new_center = true;
			}
		} 

	}

	// if (found_new_center) {
	// 	std::cout << "found new center\n";
	// }
	// else {
	// 	std::cout << "no new center... :( \n";
	// }
	return found_new_center;

}

bool ensembleMethod1(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {

	if (naiveCenterUpdate(curves, clustering, dist_func, c2c_dist)) {
		std::cout << "";
	}

	// bool new_centers = calcFSACenters(curves, clustering, 0, dist_func, c2c_dist, CenterCurveUpdateMethod::frechetMean);
	bool new_centers = cdba(curves, clustering, dist_func, c2c_dist);
	return new_centers;
}
