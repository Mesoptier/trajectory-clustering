#pragma once

// #include "../geom.h"
// #include "../IntegralFrechet/IntegralFrechet.h"
// #include <map>
// #include "util/wedge.h"
// #include "util/linear_regression.h"

// using Lines = std::vector<Line>;

// bool regression_method(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {

//     bool found_new_center = false;

//     for (auto& cluster: clustering) {
// 		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
// 			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
// 		}
// 	}

//     for (auto& cluster: clustering) {
//         auto center_curve = cluster.center_curve;
//         Lines reg_lines = Lines();
//         std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

//         for (size_t i = 1; i < center_curve.size(); ++i) {

//             Points points = Points();
//             std::vector<distance_t> weights = std::vector<distance_t>();


//             for (auto& curve_id: cluster.curve_ids) {
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

//                 WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, i - 1, 0);

//                 for (auto wp: wps) {
//                     points.push_back(wp.point);
//                     weights.push_back(wp.weight);
//                 }
//             }

//             reg_lines.push_back(
//                 linear_regression(points, weights)
//             );
//         }

//         Points new_vertices = Points();

//         new_vertices.push_back(
//             reg_lines[0].closest(center_curve[0])
//         );

//         for (int i = 1; i < reg_lines.size(); ++i) {
//             new_vertices.push_back(
//                 intersect(reg_lines[i], reg_lines[i-1])
//             );
//         }

//         new_vertices.push_back(
//             reg_lines.back().closest(center_curve.back())
//         );

//         Curve new_center_curve = Curve(new_vertices);

//         if (center_curve != new_center_curve) {
// 			// auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
//             // if (new_dist < cluster.cost) {
// 				// std::cout << "new cluster center\n";
// 				cluster.center_curve = std::move(new_center_curve);
// 				// cluster.cost = new_dist;
// 				found_new_center = true;
// 			// }
// 		}
//     }

//     if (found_new_center) {
// 		std::cout << "found new center\n";
// 	}
// 	else {
// 		std::cout << "no new center... :( \n";
// 	}

//     return false;
// }

// bool regression_method_3d(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
//        bool found_new_center = false;

//     for (auto& cluster: clustering) {
// 		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
// 			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
// 		}
// 	}

//     for (auto& cluster: clustering) {
//         auto center_curve = cluster.center_curve;
//         Lines reg_lines = Lines();
//         std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

//         // PLOTTING
//         std::fstream script;
//         script.open("reg_plot.txt", std::fstream::out | std::fstream::trunc);
//         script << "plot ";

//         std::vector<std::string> colors = {
//             "red", "blue", "green", "orange", "purple", "yellow", "grey", "cyan", "black", "brown"
//         };

//         for (size_t i = 1; i < center_curve.size(); ++i) {

//             Points points = Points();
//             std::vector<distance_t> weights = std::vector<distance_t>();
//             std::vector<distance_t> heights = std::vector<distance_t>();

//             std::fstream point_file;
//             point_file.open("points/points_" + std::to_string(i-1) + ".txt", std::fstream::out | std::fstream::trunc);    

//             for (auto& curve_id: cluster.curve_ids) {
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

//                 WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, i-1, 0);

//                 for (auto wp: wps) {

//                     point_file << wp.point.x << " " << wp.point.y << "\n\n";

//                     points.push_back(wp.point);
//                     weights.push_back(wp.weight);
//                     heights.push_back(wp.height);
//                 }
//             }

//             point_file.close();
//             script << "\"" << "points/points_" + std::to_string(i-1) + ".txt" + "\" with linespoints ls 1 lw 1 lt rgb \"" + colors[i-1] + "\", ";

//             reg_lines.push_back(
//                 linear_regression_3d(points, weights, heights)
//             );
//         }

//         for (size_t i = 0; i < center_curve.size() - 1; ++i) {
//             std::fstream segment_file;
//             segment_file.open("segments/segment_" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
//             Point source, target;
//             if (i > 0 && i < center_curve.size()-2) {
//                 source = intersect(reg_lines[i-1], reg_lines[i]);
//                 target = intersect(reg_lines[i], reg_lines[i+1]);
//             }

//             if (i == 0) {
//                 source = reg_lines[0].closest(center_curve[0]);
//                 target = intersect(reg_lines[i], reg_lines[i+1]);
//             }

//             if (i == center_curve.size()-2) {
//                 source = intersect(reg_lines[i-1], reg_lines[i]);
//                 target = reg_lines.back().closest(center_curve.back());
//             }

//             segment_file << source.x << " " << source.y << "\n" << target.x << " " << target.y;
//             segment_file.close();
//             script << "\"" << "segments/segment_" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 4 lt rgb \"" + colors[i] + "\", ";
//         }
        
//         script.close();
//         Points new_vertices = Points();

//         new_vertices.push_back(
//             reg_lines[0].closest(center_curve[0])
//         );

//         for (int i = 1; i < reg_lines.size(); ++i) {
//             new_vertices.push_back(
//                 intersect(reg_lines[i], reg_lines[i-1])
//             );
//         }

//         new_vertices.push_back(
//             reg_lines.back().closest(center_curve.back())
//         );

//         Curve new_center_curve = Curve(new_vertices);

//         if (center_curve != new_center_curve) {
// 			// auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist, dist_func);
//             // if (new_dist < cluster.cost) {
// 				// std::cout << "new cluster center\n";
// 				cluster.center_curve = std::move(new_center_curve);
// 				// cluster.cost = new_dist;
// 				found_new_center = true;
// 			// }
// 		}
//     }

//     if (found_new_center) {
// 		std::cout << "found new center\n";
// 	}
// 	else {
// 		std::cout << "no new center... :( \n";
// 	}

//     return false;
// }