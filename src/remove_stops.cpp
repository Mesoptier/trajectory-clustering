#include "remove_stops.h"


Points interpolate_between_endpoints(Points& points) {
    Point direction_vector = points.back() - points.front();
    direction_vector = direction_vector / direction_vector.dist({0, 0});
    distance_t dist = (points.back().dist(points.front())) / (points.size()-1);

    Points new_points = Points();

    for (int i = 0; i < points.size(); ++i) {
        new_points.push_back(
            points[0] + direction_vector * i
        );
    }

    return new_points;
}

Curve remove_stops(Curve curve) {

    int partition_size = 20;

    std::vector<Points> partition = std::vector<Points>();
    int next_point_index = 0;

    for (int i = 0; i < partition_size; ++i) {
        partition.push_back(Points());
        while (next_point_index < curve.size() 
        && partition.back().size() < curve.size() / partition_size) {
            partition.back().push_back(curve[next_point_index]);
            next_point_index++;
        }
    }

    std::vector<distance_t> radii = std::vector<distance_t>();

    distance_t average_radius = 0;

    for (auto& points: partition) {
        radii.push_back(
            calcMinEnclosingCircle(points).radius
        );

        average_radius += radii.back();
    }

    average_radius /= partition_size;
    distance_t threshold_ratio = 0.5;

    Curve stop_free_curve = Curve();

    for (int i = 0; i < partition.size(); ++i) {
        if (radii[i] >= average_radius * threshold_ratio) {
            for (auto point: partition[i]) {
                stop_free_curve.push_back(point);
            }
        } else {
            // std::cout << "removing stops\n";
            Points new_points = interpolate_between_endpoints(partition[i]);
            for (auto point: new_points) {
                stop_free_curve.push_back(point);
            }
        }
    }

    return stop_free_curve;
}