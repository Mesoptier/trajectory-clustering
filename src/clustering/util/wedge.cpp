#include "clustering/util/wedge.h"

#include <algorithm>
#include <iterator>

auto clustering::get_points_matched_to_segment(Points const& param_space_path,
        Curve const& curve_1, Curve const& curve_2,
        std::size_t src_index, unsigned seg_index) -> WedgePoints {
    assert(src_index < curve_1.size() - 1);
    // assert(seg_index == 0 || seg_index == 1);

    WedgePoints wedge_points;
    distance_t src_dist = curve_1.curve_length(src_index);
    distance_t tgt_dist = curve_1.curve_length(src_index + 1);

    std::vector<distance_t> x_coords(param_space_path.size());
    std::size_t p = 0;
    for (auto const& point: param_space_path) {
        x_coords[p] = point.x;
        ++p;
    }

    auto it = std::lower_bound(x_coords.begin(), x_coords.end(), src_dist);
    auto src_ind = static_cast<std::size_t>(std::distance(x_coords.begin(), it));
    it = std::lower_bound(it, x_coords.end(), tgt_dist);
    auto tgt_ind = static_cast<std::size_t>(std::distance(x_coords.begin(), it));
    if (tgt_ind == param_space_path.size())
        --tgt_ind;

    distance_t src_y = param_space_path[src_ind].y;
    distance_t tgt_y = param_space_path[tgt_ind].y;

    auto start_ind = curve_2.get_cpoint_after(src_y).ceil().getPoint();
    auto end_ind =  curve_2.get_cpoint_after(tgt_y).ceil().getPoint();

    for (auto i = start_ind; i <= end_ind; ++i) {
        distance_t weight = 0.0;
        if (i > 0) 
            weight += curve_2[i].dist(curve_2[i - 1]);
        if (i < curve_2.size() - 1)
            weight += curve_2[i].dist(curve_2[i + 1]);
        wedge_points.emplace_back(curve_2[i], seg_index, weight);
    }
    return wedge_points;
}

distance_t clustering::wedge_cost(Wedge const& wedge) {
    distance_t cost = 0.0;
    for (auto const& wp: wedge.wedge_points) {
        unsigned i = wp.matching_segment_index;
        auto dist = segPointDist(wedge.vertices[i], wedge.vertices[i + 1],
            wp.point);
        cost += dist * dist * wp.weight;
    }
    return cost;
}

Point clustering::grid_search(Wedge wedge, distance_t eps, int radius) {
    distance_t min_cost = std::numeric_limits<distance_t>::max();
    Point best = wedge.vertices[1];

    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            Point new_point{wedge.vertices[1].x + eps * i,
                wedge.vertices[1].y + eps * j};
            wedge.vertices[1] = new_point;
            if (!approx_equal(wedge.vertices[0], wedge.vertices[1])
                    && !approx_equal(wedge.vertices[1], wedge.vertices[2])) {
                auto new_cost = wedge_cost(wedge);
                if (new_cost < min_cost) {
                    min_cost = new_cost;
                    best = new_point;
                }
            }
        }
    }
    return best;
}
