#include "center_update.h"
#include <math.h>
#include <stdio.h>
#include "../IntegralFrechet/IntegralFrechet.h"
#include "../Frechet/frechet_light.h"
#include "../Frechet/frechet_matching.h"
#include "../DTW/dtw.h"
#include "../geom.h"
#include <map>
#include "util/wedge.h"
#include "util/linear_regression.h"

using Lines = std::vector<Line>;

namespace {
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

            if (approx_equal(num, 0.))
                num = 0.;

            return {num, num};
        }

        int next_index = index + 1;

        while (next_index < param_space_path.size() && param_space_path[next_index].x == param_space_path[index].x)
            ++next_index;
        

        if (isnan(param_space_path[index].y))
            std::cout << "uh-oh spaghetti-ohs\n";
        
        return {param_space_path[index].y, param_space_path[next_index - 1].y};
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

        Point new_point = Point(x_mean, y_mean);

        return new_point;
    }

    distance_t wedge_cost(Points& wedge, WedgePoints& wedge_points) {
        // std::cout << "cost start...\n";
        // std::cout << wedge[0] << ", " << wedge[1] << ", " << wedge[2] << "\n";
        // assert(wedge.size() == 3);
        
        distance_t cost = 0;

	    for (auto& wp: wedge_points) {
            
            // distance_t cost_1 = segPointDist(wedge[0], wedge[1], wp.point) 
            // * segPointDist(wedge[0], wedge[1], wp.point) 
            // * wp.weight;
            
            // distance_t cost_2 = segPointDist(wedge[1], wedge[2], wp.point)
            // * segPointDist(wedge[1], wedge[2], wp.point)
            // * wp.weight;
            
            // cost += std::min(cost_1, cost_2);

            if (wp.matching_segment_index == 0) {
                cost += segPointDist(wedge[0], wedge[1], wp.point) 
                * segPointDist(wedge[0], wedge[1], wp.point)
                * wp.weight;
            } else if (wp.matching_segment_index == 1) {
			    cost += segPointDist(wedge[1], wedge[2], wp.point) 
                * segPointDist(wedge[1], wedge[2], wp.point)
                * wp.weight;
		    } else if (wp.matching_segment_index == 2) {
                cost += segPointDist(wedge[2], wedge[3], wp.point) 
                * segPointDist(wedge[2], wedge[3], wp.point)
                * wp.weight;
            }
	    }
        // std::cout << "cost end...\n";
	    return cost;
    }

    Point grid_search(Points wedge, WedgePoints wedge_points, distance_t eps = 0.125, int radius = 20) {
    
        Points candidate_wedge = {wedge[0], wedge[1], wedge[2]};
        distance_t min_cost = std::numeric_limits<distance_t>::max();
        Point best = wedge[1];


        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                Point new_point = {wedge[1].x + eps * i, wedge[1].y + eps * j};
                candidate_wedge[1] = new_point;
                if (!approx_equal(wedge[0], wedge[1]) && !approx_equal(wedge[1], wedge[2])) {
                    distance_t new_cost = wedge_cost(candidate_wedge, wedge_points);
                    if (new_cost < min_cost) {
                        min_cost = new_cost;
                        best = new_point;
                    }
                }
            }
        }

        return best;
    }


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
}

Curve dba_update(Curves const& curves, Cluster const& cluster) {
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

    new_center_curve.push_back(cluster.center_curve[0]);

	for (int i = 1; i < cluster.center_curve.size(); ++i) {
		Points points_to_average = Points();

		for (auto& matching: matchings) {
			// points_to_average.push_back(matching[i][0]);
			for (auto& point: matching[i]) {
				points_to_average.push_back(point);
			}
		}

		Point new_point = mean_of_points(points_to_average);

		if (new_center_curve.size() == 0 || !approx_equal(new_point, new_center_curve.back())) {
			new_center_curve.push_back(new_point);
		}
    }

    //Fix last point for pigeon data
    // new_center_curve.push_back(cluster.center_curve.back());

    return new_center_curve;
}

Curve cdba_update(Curves const& curves, Cluster const& cluster) {
        
    const auto& center_curve = cluster.center_curve;
	Curve new_center_curve;
	std::vector<std::vector<Points>> matchings = std::vector<std::vector<Points>>();

	for (auto& curve_id: cluster.curve_ids) {
		Curve curve = curves[curve_id];
		Points param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 250, nullptr)
		.compute_matching()
		.matching;

		std::vector<Points> matching = std::vector<Points>();

		for (int i = 0; i < center_curve.size(); ++i) {
			matching.push_back({});
			distance_t dist = center_curve.curve_length(i);
			std::pair<distance_t, distance_t> y_range = get_y_range(param_space_path, dist);

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
		}

		matchings.push_back(matching);
	}

    new_center_curve.push_back(cluster.center_curve[0]);

	for (int i = 1; i < cluster.center_curve.size() ; ++i) {
		Points points_to_average = Points();

		for (auto& matching: matchings) {
			for (auto& point: matching[i]) {
				points_to_average.push_back(point);
			}
		}

		Point new_point = mean_of_points(points_to_average);

		if (new_center_curve.size() == 0 || !approx_equal(new_point, new_center_curve.back())) {
			new_center_curve.push_back(new_point);
		}
	}

    // Fix last point for pigeon data
    // new_center_curve.push_back(cluster.center_curve.back());

    return new_center_curve;
}

Curve cdba_update(Curves const& curves, Cluster const& cluster, int res) {
        
    const auto& center_curve = cluster.center_curve;
	Curve new_center_curve;
	std::vector<std::vector<Points>> matchings = std::vector<std::vector<Points>>();

	for (auto& curve_id: cluster.curve_ids) {
		Curve curve = curves[curve_id];
		Points param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 5, nullptr)
		.compute_matching()
		.matching;

		std::vector<Points> matching = std::vector<Points>();

		for (int i = 0; i < center_curve.size(); ++i) {
			matching.push_back({});
			distance_t dist = center_curve.curve_length(i);
			std::pair<distance_t, distance_t> y_range = get_y_range(param_space_path, dist);

			if (approx_equal(y_range.first, y_range.second)) {
				matching.back().push_back(curve.interpolate_at(curve.get_cpoint_after(y_range.first)));
			}
			else {
				distance_t length = y_range.second - y_range.first;
				distance_t ratio_to_curve_length = length / curve.curve_length();
				int number_of_samples = res * curve.size() * ratio_to_curve_length;

				if (number_of_samples == 0)
					number_of_samples = 1;

				for (int j = 0; j < number_of_samples; ++j) {
					distance_t dist_along_curve = y_range.first + j * length / number_of_samples;
					matching.back().push_back(curve.interpolate_at(curve.get_cpoint_after(dist_along_curve)));
				}
			}
		}

		matchings.push_back(matching);
	}

    new_center_curve.push_back(cluster.center_curve[0]);

	for (int i = 1; i < cluster.center_curve.size(); ++i) {
		Points points_to_average = Points();

		for (auto& matching: matchings) {
			// points_to_average.push_back(matching[i][0]);
			for (auto& point: matching[i]) {
				points_to_average.push_back(point);
			}
		}


		Point new_point = mean_of_points(points_to_average);

		if (new_center_curve.size() == 0 || !approx_equal(new_point, new_center_curve.back())) {
			new_center_curve.push_back(new_point);
		}
	}

    //Fix last point for pigeon data
    // new_center_curve.push_back(cluster.center_curve.back());

    return new_center_curve;
}

Curve fsa_update(Curves const& curves, Cluster const& cluster) {

	std::vector<Points> matchings;
		auto const& center_curve = cluster.center_curve;
		Curve new_center_curve;

		for (auto curve_id: cluster.curve_ids) {
			auto const& curve = curves[curve_id];
			auto matching = calcMatching(cluster.center_curve, curve);
			matchings.push_back(std::move(matching));
		}

		for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
			Points matching_points;
			for (auto const& matching: matchings) {
				matching_points.push_back(matching[point_id]);
			}
			auto min_enclosing_circle = calcMinEnclosingCircle(matching_points);
            Point new_point = min_enclosing_circle.center;
            if (new_center_curve.get_points().empty() || !approx_equal(new_point, new_center_curve.get_points().back())) {
				new_center_curve.push_back(new_point);
			}
        }

    return new_center_curve;
}

Curve regression_update(Curves const& curves, Cluster const& cluster) {
	auto center_curve = cluster.center_curve;
    Lines reg_lines = Lines();
    std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

    for (size_t i = 1; i < center_curve.size(); ++i) {

        Points points = Points();
        std::vector<distance_t> weights = std::vector<distance_t>();


        for (auto& curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];
            Points param_space_path;

            if (matching_paths.find(curve_id) == matching_paths.end()) {
                param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 5, nullptr)
                .compute_matching()
                .matching;
                matching_paths.emplace(curve_id, param_space_path);
            } else {
                param_space_path = matching_paths.at(curve_id);
            }

            WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, i - 1, 0);

            for (auto wp: wps) {
                points.push_back(wp.point);
                weights.push_back(wp.weight);
            }
        }

        reg_lines.push_back(
            linear_regression(points, weights)
        );
    }

    Points new_vertices = Points();

    new_vertices.push_back(
        reg_lines[0].closest(center_curve[0])
    );

    for (int i = 1; i < reg_lines.size(); ++i) {
        new_vertices.push_back(
            intersect(reg_lines[i], reg_lines[i-1])
        );
    }

    new_vertices.push_back(
        reg_lines.back().closest(center_curve.back())
    );

    Curve new_center_curve = Curve(new_vertices);

	return new_center_curve;
}

Curve regression_update_3d(Curves const& curves, Cluster const& cluster) {
    std::cout << "computing new center curve...\n";
	auto center_curve = cluster.center_curve;
    Lines reg_lines = Lines();
    std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

    // PLOTTING
    std::fstream script;
    script.open("reg_plot.txt", std::fstream::out | std::fstream::trunc);
    script << "plot ";

    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "yellow", "grey", "cyan", "black", "brown"
    };

    for (size_t i = 1; i < center_curve.size(); ++i) {

        Points points = Points();
        std::vector<distance_t> weights = std::vector<distance_t>();
        std::vector<distance_t> heights = std::vector<distance_t>();

        std::fstream point_file;
        point_file.open("points/points_" + std::to_string(i-1) + ".txt", std::fstream::out | std::fstream::trunc);    

        for (auto& curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];
            Points param_space_path;

            if (matching_paths.find(curve_id) == matching_paths.end()) {
                param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 5, nullptr)
                .compute_matching()
                 .matching;
                matching_paths.emplace(curve_id, param_space_path);
            } else {
                param_space_path = matching_paths.at(curve_id);
            }

            WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, i-1, 0);

            for (auto wp: wps) {

                point_file << wp.point.x << " " << wp.point.y << "\n\n";

                points.push_back(wp.point);
                weights.push_back(wp.weight);
                heights.push_back(wp.height);
            }
        }

        point_file.close();
        script << "\"" << "points/points_" + std::to_string(i-1) + ".txt" + "\" with linespoints ls 1 lw 1 lt rgb \"" + colors[i-1] + "\", ";

        reg_lines.push_back(
            linear_regression_3d(points, weights, heights)
        );
    }

    for (size_t i = 0; i < center_curve.size() - 1; ++i) {
        std::fstream segment_file;
        segment_file.open("segments/segment_" + std::to_string(i) + ".txt", std::fstream::out | std::fstream::trunc);
        Point source, target;
        if (i > 0 && i < center_curve.size()-2) {
            source = intersect(reg_lines[i-1], reg_lines[i]);
            target = intersect(reg_lines[i], reg_lines[i+1]);
        }

        if (i == 0) {
            source = reg_lines[0].closest(center_curve[0]);
            target = intersect(reg_lines[i], reg_lines[i+1]);
        }

        if (i == center_curve.size()-2) {
            source = intersect(reg_lines[i-1], reg_lines[i]);
            target = reg_lines.back().closest(center_curve.back());
        }

        segment_file << source.x << " " << source.y << "\n" << target.x << " " << target.y;
        segment_file.close();
        script << "\"" << "segments/segment_" + std::to_string(i) + ".txt" + "\" with linespoints ls 1 lw 4 lt rgb \"" + colors[i] + "\", ";
    }
        
    script.close();
    Points new_vertices = Points();

    new_vertices.push_back(
        reg_lines[0].closest(center_curve[0])
    );

    for (int i = 1; i < reg_lines.size(); ++i) {
        new_vertices.push_back(
            intersect(reg_lines[i], reg_lines[i-1])
        );
    }

    new_vertices.push_back(
        reg_lines.back().closest(center_curve.back())
    );

    Curve new_center_curve = Curve(new_vertices);
	return new_center_curve;
}

Curve wedge_update(Curves const& curves, Cluster const& cluster) {

	// distance_t min_edge_length = std::numeric_limits<distance_t>::max();
	// distance_t max_edge_length = std::numeric_limits<distance_t>::min();	
    // for (auto& curve: curves) {
    //     for (size_t i = 1; i < curve.size(); ++i) {
    //         distance_t dist = curve.curve_length(i-1, i);
    //         if (dist < min_edge_length) {
    //             min_edge_length = dist;
    //         }
    //         if (dist > max_edge_length) {
    //             max_edge_length = dist;
    //         }
    //     }
    // }
	

    const auto& center_curve = cluster.center_curve;
    Curve new_center_curve;
    Wedges wedges = Wedges();
    new_center_curve.push_back(center_curve[0]);

    std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

    for (size_t i = 0; i < center_curve.size() - 2; ++i) {
        Points vertices = new_center_curve.size() == 0 
        ? Points({center_curve[i], center_curve[i+1], center_curve[i+2]})
        : Points({new_center_curve.back(), center_curve[i+1], center_curve[i+2]});

        wedges.push_back(Wedge(
            vertices, WedgePoints()
        ));

        for (auto curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];

            Points param_space_path;

            if (matching_paths.find(curve_id) == matching_paths.end()) {
                param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 250, nullptr)
                .compute_matching()
                .matching;
                matching_paths.emplace(curve_id, param_space_path);
            } else {
                param_space_path = matching_paths.at(curve_id);
            }

            WedgePoints seg_1_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i, 0);
            WedgePoints seg_2_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i + 1, 1);

            if (wedges.back().wedge_points.empty()) {
                wedges.back().wedge_points = seg_1_points;
            } else {
                wedges.back().wedge_points.insert(
                    wedges.back().wedge_points.end(), seg_1_points.begin(), seg_1_points.end()
                );
            }
            wedges.back().wedge_points.insert(
                wedges.back().wedge_points.end(), seg_2_points.begin(), seg_2_points.end()
            );
        }

        Point new_point = grid_search(wedges.back().vertices, wedges.back().wedge_points/*, 5, 10*/);

        if (!approx_equal(new_point, new_center_curve.back())) {
            new_center_curve.push_back(new_point);
        }
    }

    Points last_points = {center_curve.get_points()[center_curve.size() - 2], center_curve.get_points()[center_curve.size() - 1], {0, 0}};
    Wedge last_wedge = Wedge(
        last_points,
        WedgePoints()
    );
    for (auto curve_id: cluster.curve_ids) {
        Curve curve = curves[curve_id];
        Points param_space_path = matching_paths.at(curve_id);
        WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, center_curve.size() - 2, 0);
        last_wedge.wedge_points.insert(last_wedge.wedge_points.end(), wps.begin(), wps.end());
    }

    Point new_point = grid_search(last_wedge.vertices, last_wedge.wedge_points);

    if (new_center_curve.empty() || !approx_equal(new_point, new_center_curve.back())) {
        new_center_curve.push_back(new_point);
    }

    // new_center_curve.push_back(center_curve.back());

    return new_center_curve;
}

Curve wedge_update_param_args(Curves const& curves, Cluster const& cluster, distance_t eps, int radius) {

    const auto& center_curve = cluster.center_curve;
    Curve new_center_curve;
    Wedges wedges = Wedges();
    new_center_curve.push_back(center_curve[0]);

    std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

    for (size_t i = 0; i < center_curve.size() - 2; ++i) {
        Points vertices = new_center_curve.size() == 0 
        ? Points({center_curve[i], center_curve[i+1], center_curve[i+2]})
        : Points({new_center_curve.back(), center_curve[i+1], center_curve[i+2]});

        wedges.push_back(Wedge(
            vertices, WedgePoints()
        ));

        for (auto curve_id: cluster.curve_ids) {
            Curve curve = curves[curve_id];

            Points param_space_path;

            if (matching_paths.find(curve_id) == matching_paths.end()) {
                param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 5, nullptr)
                .compute_matching()
                .matching;
                matching_paths.emplace(curve_id, param_space_path);
            } else {
                param_space_path = matching_paths.at(curve_id);
            }

            WedgePoints seg_1_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i, 0);
            WedgePoints seg_2_points = get_points_matched_to_segment(param_space_path, center_curve, curve, i + 1, 1);

            if (wedges.back().wedge_points.empty()) {
                wedges.back().wedge_points = seg_1_points;
            } else {
                wedges.back().wedge_points.insert(
                    wedges.back().wedge_points.end(), seg_1_points.begin(), seg_1_points.end()
                );
            }
            wedges.back().wedge_points.insert(
                wedges.back().wedge_points.end(), seg_2_points.begin(), seg_2_points.end()
            );
        }

        Point new_point = grid_search(wedges.back().vertices, wedges.back().wedge_points, eps, radius);

        if (!approx_equal(new_point, new_center_curve.back())) {
            new_center_curve.push_back(new_point);
        }
    }

    Points last_points = {center_curve.get_points()[center_curve.size() - 2], center_curve.get_points()[center_curve.size() - 1], {0, 0}};
    Wedge last_wedge = Wedge(
        last_points,
        WedgePoints()
    );
    for (auto curve_id: cluster.curve_ids) {
        Curve curve = curves[curve_id];
        Points param_space_path = matching_paths.at(curve_id);
        WedgePoints wps = get_points_matched_to_segment(param_space_path, center_curve, curve, center_curve.size() - 2, 0);
        last_wedge.wedge_points.insert(last_wedge.wedge_points.end(), wps.begin(), wps.end());
    }

    Point new_point = grid_search(last_wedge.vertices, last_wedge.wedge_points);

    if (new_center_curve.empty() || !approx_equal(new_point, new_center_curve.back())) {
        new_center_curve.push_back(new_point);
    }

    new_center_curve.back() == center_curve.back();

    return new_center_curve;
}

Curve redistribute_points_update(Curves const& curves, Cluster const& cluster) {
    std::vector<Points> matchings;
	auto const& center_curve = cluster.center_curve;
	Curve new_center_curve;

	for (auto curve_id: cluster.curve_ids) {
		auto const& curve = curves[curve_id]; 
        matchings.push_back(
            get_images_of_points(
                IntegralFrechet(cluster.center_curve, curve, ParamMetric::LInfinity_NoShortcuts, 5, nullptr)
                .compute_matching()
                .matching,
                get_redistributed_distances(cluster.center_curve),
                curve
            )
        );
    }


    for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
        Points matching_points;
        for (auto const& matching: matchings) {
            matching_points.push_back(matching[point_id]);
        }

        Point new_point = mean_of_points(matching_points);

        if (new_center_curve.get_points().empty() || !approx_equal(new_point, new_center_curve.get_points().back())) {
            new_center_curve.push_back(new_point);
        }
    }

    return new_center_curve;
}

