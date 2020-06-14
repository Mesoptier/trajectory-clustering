#pragma once
// #include "../geom.h"
// #include "../IntegralFrechet/IntegralFrechet.h"
// #include "util/wedge.h"
// #include <map>



// Point mean_of_points(Points points) {

// 	// std::cout << "\n";
// 	// for (auto& point: points)
// 	// 	std::cout << point << ", ";
// 	// std::cout << "\n";

// 	distance_t x_mean = 0;
// 	distance_t y_mean = 0;

// 	for (auto p: points) {
// 		x_mean += p.x;
// 		y_mean += p.y;
// 	}

// 	x_mean /= points.size();
// 	y_mean /= points.size();

// 	Point new_point = Point(x_mean, y_mean);

// 	return new_point;
// }

// distance_t wedge_cost(Points& wedge, WedgePoints& wedge_points) {
// 	// std::cout << "cost start...\n";
//     // std::cout << wedge[0] << ", " << wedge[1] << ", " << wedge[2] << "\n";
// 	// assert(wedge.size() == 3);
	
// 	distance_t cost = 0;

// 	for (auto& wp: wedge_points) {
		
//         // distance_t cost_1 = segPointDist(wedge[0], wedge[1], wp.point) 
//         // * segPointDist(wedge[0], wedge[1], wp.point) 
//         // * wp.weight;
        
//         // distance_t cost_2 = segPointDist(wedge[1], wedge[2], wp.point)
//         // * segPointDist(wedge[1], wedge[2], wp.point)
//         // * wp.weight;
        
//         // cost += std::min(cost_1, cost_2);

//         if (wp.matching_segment_index == 0) {
// 			cost += segPointDist(wedge[0], wedge[1], wp.point) 
//             * segPointDist(wedge[0], wedge[1], wp.point)
//             * wp.weight;
// 		} else if (wp.matching_segment_index == 1) {
// 			cost += segPointDist(wedge[1], wedge[2], wp.point) 
//             * segPointDist(wedge[1], wedge[2], wp.point)
//             * wp.weight;
// 		} else if (wp.matching_segment_index == 2) {
//             cost += segPointDist(wedge[2], wedge[3], wp.point) 
//             * segPointDist(wedge[2], wedge[3], wp.point)
//             * wp.weight;
//         }
// 	}
//     // std::cout << "cost end...\n";
// 	return cost;
// }


// Point grid_search(Points wedge, WedgePoints wedge_points, distance_t eps = 0.25) {
    
//     Points candidate_wedge = {wedge[0], wedge[1], wedge[2]};
//     distance_t min_cost = std::numeric_limits<distance_t>::max();
//     Point best = wedge[1];


//     for (int i = -10; i <= 10; ++i) {
//         for (int j = -10; j <= 10; ++j) {
//             Point new_point = {wedge[1].x + eps * i, wedge[1].y + eps * j};
//             candidate_wedge[1] = new_point;
//             if (!approx_equal(wedge[0], wedge[1]) && !approx_equal(wedge[1], wedge[2])) {
//                 distance_t new_cost = wedge_cost(candidate_wedge, wedge_points);
//                 if (new_cost < min_cost) {
//                     min_cost = new_cost;
//                     best = new_point;
//                 }
//             }
//         }
//     }

//     return best;
// }

// std::pair<Point, Point> grid_search_2 (Points wedge, WedgePoints wedge_points, distance_t eps = 0.125) {
//     Points candidate_wedge = {wedge[0], wedge[1], wedge[2], wedge[3]};
//     distance_t min_cost = std::numeric_limits<distance_t>::max();
//     Point best_1 = wedge[1];
//     Point best_2 = wedge[2];


//     for (int i = -5; i <= 5; ++i) {
//         for (int j = -5; j <= 5; ++j) {
//             for (int k = -5; k <= 5; ++k) {
//                 for (int l = -5; l <= 5; ++l) {
//                     Point new_point_1 = {wedge[1].x + eps * i, wedge[1].y + eps * j};
//                     Point new_point_2 = {wedge[2].x + eps * k, wedge[2].y + eps * l};
//                     candidate_wedge[1] = new_point_1;
//                     candidate_wedge[2] = new_point_2;
//                     if (!approx_equal(wedge[0], wedge[1]) && 
//                     !approx_equal(wedge[1], wedge[2]) && 
//                     !approx_equal(wedge[2], wedge[3])) {
//                         distance_t new_cost = wedge_cost(candidate_wedge, wedge_points);
//                         if (new_cost < min_cost) {
//                             min_cost = new_cost;
//                             best_1 = new_point_1;
//                             best_2 = new_point_2;
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return {best_1, best_2};
// }
// bool wedge_method_2 (Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
//     bool found_new_center = false;

//     for (auto& cluster: clustering) {
// 		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
// 			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
// 		}
// 	}

//     for (auto& cluster: clustering) {
// 		const auto& center_curve = cluster.center_curve;
// 		Curve new_center_curve;
// 		Wedges wedges = Wedges();
//         new_center_curve.push_back(center_curve[0]);

//         std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

// 		for (size_t i = 0; i < center_curve.size() - 3; ++i) {
//             Points vertices = {center_curve[i], center_curve[i+1], center_curve[i+2], center_curve[i+3]};

//             wedges.push_back(Wedge(
//                 vertices, WedgePoints()
//             ));

//             for (auto curve_id: cluster.curve_ids) {
//                 Curve curve = curves[curve_id];

//                 Points param_space_path;

//                 if (matching_paths.find(curve_id) == matching_paths.end()) {
//                     param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 1, nullptr)
//                     .compute_matching()
//                     .matching;
//                     matching_paths.emplace(curve_id, param_space_path);
//                 } else {
//                     param_space_path = matching_paths.at(curve_id);
//                 }

//                 WedgePoints seg_1_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i, 0);
//                 WedgePoints seg_2_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i + 1, 1);
//                 WedgePoints seg_3_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i + 1, 2);

//                 if (wedges.back().wedge_points.empty()) {
//                     wedges.back().wedge_points = seg_1_points;
//                 } else {
//                     wedges.back().wedge_points.insert(
//                         wedges.back().wedge_points.end(), seg_1_points.begin(), seg_1_points.end()
//                     );
//                 }
//                 wedges.back().wedge_points.insert(
//                     wedges.back().wedge_points.end(), seg_2_points.begin(), seg_2_points.end()
//                 );
//                 wedges.back().wedge_points.insert(
//                     wedges.back().wedge_points.end(), seg_3_points.begin(), seg_3_points.end()
//                 );
//             }

//             std::pair<Point, Point> new_points = grid_search_2(wedges.back().vertices, wedges.back().wedge_points);

//             if (!approx_equal(new_points.first, new_center_curve.back())) {
//                 new_center_curve.push_back(new_points.first);
//             }

//             if (!approx_equal(new_points.second, new_center_curve.back())) {
//                 new_center_curve.push_back(new_points.second);
//             }
//         }

//         Points last_points = {center_curve.get_points()[center_curve.size() - 2], center_curve.get_points()[center_curve.size() - 1], {0, 0}};
//         Wedge last_wedge = Wedge(
//             last_points,
//             WedgePoints()
//         );
//         for (auto curve_id: cluster.curve_ids) {
//             Curve curve = curves[curve_id];
//             Points param_space_path = matching_paths.at(curve_id);
//             WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, center_curve.size() - 2, 0);
//             last_wedge.wedge_points.insert(last_wedge.wedge_points.end(), wps.begin(), wps.end());
//         }

//         Point new_point = grid_search(last_wedge.vertices, last_wedge.wedge_points);

//         // Points last_points = Points();
//         // for (auto curve_id: cluster.curve_ids) {
//         //     last_points.push_back(curves[curve_id].get_points().back());
//         // }

//         // Point new_point = mean_of_points(last_points);

//         if (new_center_curve.empty() || !approx_equal(new_point, new_center_curve.back())) {
//             new_center_curve.push_back(new_point);
//         }


// 		if (center_curve != new_center_curve) {
// 			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
//             if (new_dist < cluster.cost) {
// 				// std::cout << "new cluster center\n";
// 				cluster.center_curve = std::move(new_center_curve);
// 				cluster.cost = new_dist;
// 				found_new_center = true;
// 			}
// 		} 
// 	}

// 	if (found_new_center) {
// 		std::cout << "found new center\n";
// 	}
// 	else {
// 		std::cout << "no new center... :( \n";
// 	}
// 	return found_new_center;
// }

// bool wedge_method(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
//     bool found_new_center = false;

//     for (auto& cluster: clustering) {
// 		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
// 			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
// 		}
// 	}

//     for (auto& cluster: clustering) {
// 		const auto& center_curve = cluster.center_curve;
// 		Curve new_center_curve;
// 		Wedges wedges = Wedges();
//         new_center_curve.push_back(center_curve[0]);

//         std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

// 		for (size_t i = 0; i < center_curve.size() - 2; ++i) {
//             Points vertices = new_center_curve.size() == 0 
//             ? Points({center_curve[i], center_curve[i+1], center_curve[i+2]})
//             : Points({new_center_curve.back(), center_curve[i+1], center_curve[i+2]});

//             wedges.push_back(Wedge(
//                 vertices, WedgePoints()
//             ));

//             for (auto curve_id: cluster.curve_ids) {
//                 Curve curve = curves[curve_id];

//                 Points param_space_path;

//                 if (matching_paths.find(curve_id) == matching_paths.end()) {
//                     param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 1, nullptr)
//                     .compute_matching()
//                     .matching;
//                     matching_paths.emplace(curve_id, param_space_path);
//                 } else {
//                     param_space_path = matching_paths.at(curve_id);
//                 }

//                 WedgePoints seg_1_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i, 0);
//                 WedgePoints seg_2_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i + 1, 1);

//                 if (wedges.back().wedge_points.empty()) {
//                     wedges.back().wedge_points = seg_1_points;
//                 } else {
//                     wedges.back().wedge_points.insert(
//                         wedges.back().wedge_points.end(), seg_1_points.begin(), seg_1_points.end()
//                     );
//                 }
//                 wedges.back().wedge_points.insert(
//                     wedges.back().wedge_points.end(), seg_2_points.begin(), seg_2_points.end()
//                 );
//             }

//             Point new_point = grid_search(wedges.back().vertices, wedges.back().wedge_points);

//             if (!approx_equal(new_point, new_center_curve.back())) {
//                 new_center_curve.push_back(new_point);
//             }
//         }

//         Points last_points = {center_curve.get_points()[center_curve.size() - 2], center_curve.get_points()[center_curve.size() - 1], {0, 0}};
//         Wedge last_wedge = Wedge(
//             last_points,
//             WedgePoints()
//         );
//         for (auto curve_id: cluster.curve_ids) {
//             Curve curve = curves[curve_id];
//             Points param_space_path = matching_paths.at(curve_id);
//             WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, center_curve.size() - 2, 0);
//             last_wedge.wedge_points.insert(last_wedge.wedge_points.end(), wps.begin(), wps.end());
//         }

//         Point new_point = grid_search(last_wedge.vertices, last_wedge.wedge_points);

//         // Points last_points = Points();
//         // for (auto curve_id: cluster.curve_ids) {
//         //     last_points.push_back(curves[curve_id].get_points().back());
//         // }

//         // Point new_point = mean_of_points(last_points);

//         if (new_center_curve.empty() || !approx_equal(new_point, new_center_curve.back())) {
//             new_center_curve.push_back(new_point);
//         }


// 		if (center_curve != new_center_curve) {
// 			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
//             if (new_dist < cluster.cost) {
// 				// std::cout << "new cluster center\n";
// 				cluster.center_curve = std::move(new_center_curve);
// 				cluster.cost = new_dist;
// 				found_new_center = true;
// 			}
// 		} 
// 	}

// 	if (found_new_center) {
// 		std::cout << "found new center\n";
// 	}
// 	else {
// 		std::cout << "no new center... :( \n";
// 	}
// 	return found_new_center;
// }
