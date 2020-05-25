#pragma once

#include "../../Curve.h"

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