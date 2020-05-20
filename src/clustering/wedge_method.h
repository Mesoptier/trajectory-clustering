#pragma once
#include "../geom.h"
#include "../IntegralFrechet/IntegralFrechet.h"
#include <map>

struct WedgePoint {

    Point point;
    size_t matching_segment_index;
    distance_t weight;

    WedgePoint (Point p, size_t i, distance_t w) 
    : point(p), matching_segment_index(i), weight(w) {}
     
};

using WedgePoints = std::vector<WedgePoint>;


struct Wedge {
    Points vertices;
    WedgePoints wedge_points;
    Wedge(Points vs, WedgePoints wps) 
    : vertices(vs), wedge_points(wps) {}
};

using Wedges = std::vector<Wedge>;



Point mean_of_points(Points points) {

	// std::cout << "\n";
	// for (auto& point: points)
	// 	std::cout << point << ", ";
	// std::cout << "\n";

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


WedgePoints get_points_matched_to_segment(Points& param_space_path, const Curve& curve_1, const Curve& curve_2, 
size_t src_index, size_t seg_index) {
	
	assert(src_index < curve_1.size() - 1);
    assert(seg_index == 0 || seg_index == 1);

	WedgePoints wedge_points = WedgePoints();

	distance_t src_dist = curve_1.curve_length(src_index);
	distance_t tgt_dist = curve_1.curve_length(src_index + 1);

	std::vector<distance_t> x_coords = std::vector<distance_t>();
	for (auto& point: param_space_path) {
		x_coords.push_back(point.x);
	}

	int src_ind = 0;
	int tgt_ind = 0;

	auto it = std::lower_bound(x_coords.begin(), x_coords.end(), src_dist);
	src_ind = static_cast<std::size_t>(std::distance(x_coords.begin(), it));
	it = std::lower_bound(it, x_coords.end(), tgt_dist);
	tgt_ind = static_cast<std::size_t>(std::distance(x_coords.begin(), it));

	distance_t src_y = param_space_path[src_ind].y;
	distance_t tgt_y = param_space_path[tgt_ind].y;

	size_t start_ind = curve_2.get_cpoint_after(src_y).ceil().getPoint();
	size_t end_ind =  curve_2.get_cpoint_after(tgt_y).ceil().getPoint();

	for (size_t i = start_ind; i <= end_ind; ++i) {
        distance_t weight = 0;
        if (i > 0) 
            weight += curve_2[i].dist(curve_2[i-1]);
        if (i < curve_2.size() - 1)
            weight += curve_2[i].dist(curve_2[i+1]);
        wedge_points.push_back(WedgePoint(
            curve_2[i], seg_index, weight
        ));
    }

	return wedge_points;
}

distance_t wedge_cost(Points& wedge, WedgePoints& wedge_points) {
	// std::cout << "cost start...\n";
    // std::cout << wedge[0] << ", " << wedge[1] << ", " << wedge[2] << "\n";
	assert(wedge.size() == 3);
	
	distance_t cost = 0;

	for (auto& wp: wedge_points) {
		if (wp.matching_segment_index == 0) {
			cost += segPointDist(wedge[0], wedge[1], wp.point) 
            * segPointDist(wedge[0], wedge[1], wp.point)
            * wp.weight;
		} else {
			cost += segPointDist(wedge[1], wedge[2], wp.point) 
            * segPointDist(wedge[1], wedge[2], wp.point)
            * wp.weight;
		}
	}
    // std::cout << "cost end...\n";
	return cost;
}

Point grid_search(Points wedge, WedgePoints wedge_points, distance_t eps = 0.25) {
    
    Points candidate_wedge = {wedge[0], wedge[1], wedge[2]};
    distance_t min_cost = std::numeric_limits<distance_t>::max();
    Point best = wedge[1];


    for (int i = -10; i <= 10; ++i) {
        for (int j = -10; j <= 10; ++j) {
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

bool wedge_method(Curves const& curves, Clustering& clustering, distance_t(*dist_func)(Curve, Curve), C2CDist c2c_dist) {
    std::cout << "hello from wedge method...\n";
    bool found_new_center = false;

    for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist, dist_func);
		}
	}

    for (auto& cluster: clustering) {
		const auto& center_curve = cluster.center_curve;
		Curve new_center_curve;
		Wedges wedges = Wedges();
        new_center_curve.push_back(center_curve[0]);

        std::map<CurveID, Points> matching_paths = std::map<CurveID, Points>();

		for (size_t i = 0; i < center_curve.size() - 2; ++i) {
            Points vertices = {center_curve[i], center_curve[i+1], center_curve[i+2]};

            wedges.push_back(Wedge(
                vertices, WedgePoints()
            ));

            for (auto curve_id: cluster.curve_ids) {
                Curve curve = curves[curve_id];

                Points param_space_path;

                if (matching_paths.find(curve_id) == matching_paths.end()) {
                    param_space_path = IntegralFrechet(center_curve, curve, ParamMetric::L1, 1, nullptr)
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

            Point new_point = grid_search(wedges.back().vertices, wedges.back().wedge_points);

            if (!approx_equal(new_point, new_center_curve.back())) {
                new_center_curve.push_back(new_point);
            }
        }

        Points last_points = Points();
        for (auto curve_id: cluster.curve_ids) {
            last_points.push_back(curves[curve_id].get_points().back());
        }

        Point new_point = mean_of_points(last_points);

        if (new_center_curve.empty() || !approx_equal(new_point, new_center_curve.back())) {
            new_center_curve.push_back(new_point);
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