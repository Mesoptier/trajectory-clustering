#include "utils/remove_stops.h"

namespace {
    /**
     * \brief Given a set of points, return a set of points with the same start
     * and end that linearly interpolates between them.
     * \param points The original points
     * \return Points on a line, uniformly spread, with the sam start and end as
     * the input.
     */
    Points interpolate_between_endpoints(Points const& points) {
        Point direction_vector = points.back() - points.front();
        direction_vector = direction_vector / direction_vector.dist({0, 0});
        Points new_points;
        new_points.reserve(points.size());
        for (std::size_t i = 0; i < points.size(); ++i)
            new_points.emplace_back(points[0] + direction_vector * i);
        return new_points;
    }
}

Curve io::remove_stops(Curve const& curve) {
    std::size_t partition_size = 20;
    std::size_t next_point_index = 0;
    std::vector<Points> partition(partition_size);

    for (std::size_t i = 0; i < partition_size; ++i) {
        while (next_point_index < curve.size() 
                && partition[i].size() < curve.size() / partition_size) {
            partition[i].push_back(curve[next_point_index]);
            next_point_index++;
        }
    }

    std::vector<distance_t> radii;
    radii.reserve(partition.size());
    distance_t average_radius = 0.0;
    for (auto const& points: partition) {
        radii.emplace_back(calcMinEnclosingCircle(points).radius);
        average_radius += radii.back();
    }

    average_radius /= partition_size;
    distance_t threshold_ratio = 0.5;
    Curve stop_free_curve;
    for (std::size_t i = 0; i < partition.size(); ++i) {
        if (radii[i] >= average_radius * threshold_ratio) {
            for (auto const& point: partition[i])
                stop_free_curve.push_back(point);
        } else {
            Points new_points = interpolate_between_endpoints(partition[i]);
            for (auto const& point: new_points)
                stop_free_curve.push_back(point);
        }
    }
    return stop_free_curve;
}
