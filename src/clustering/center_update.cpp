#include "clustering/center_update.h"

#include <algorithm>
#include <map>
#include "clustering/util/wedge.h"
#include "DTW/dtw.h"
#include "Frechet/frechet_light.h"
#include "Frechet/frechet_matching.h"
#include "IntegralFrechet/IntegralFrechet.h"

namespace {
    /**
     * \brief Get the range along the second curve matched to the point on the
     * first curve.
     * \param param_space_path The matching.
     * \param distance The distance along the first curve, indicating the point
     * of interest.
     * \return The range of distances along the second curve that correspond to
     * the range of points on the second curves all matched to the point at
     * `distance' on the first curve.
     */
    std::pair<distance_t, distance_t> get_y_range(
            Points const& param_space_path, distance_t distance) {
        std::vector<distance_t> x_coords;

        for (auto const& point: param_space_path)
            x_coords.emplace_back(point.x);

        if (distance > param_space_path.back().x)
            distance = param_space_path.back().x;

        auto it = std::lower_bound(x_coords.begin(), x_coords.end(), distance);
        auto index = static_cast<std::size_t>(
            std::distance(x_coords.begin(), it));

        if (x_coords[index] > distance && index > 0)
            index--;
        if (index == x_coords.size() - 1)
            return {param_space_path.back().y, param_space_path.back().y};

        auto line = Line::fromTwoPoints(param_space_path[index],
            param_space_path[index + 1]);

        if (!line.isVertical()) {
            distance_t num = line.getY(distance);
            if (approx_zero(num))
                num = 0.0;
            return {num, num};
        }

        auto next_index = index + 1;
        while (next_index < param_space_path.size() &&
                param_space_path[next_index].x == param_space_path[index].x)
            ++next_index;
        return {param_space_path[index].y, param_space_path[next_index - 1].y};
    }

    /**
     * \brief Take the arithmetic mean of x- and y-coordinates for points and
     * return the resulting point.
     * \param points The points to take the mean of.
     * \return The mean point.
     */
    Point mean_of_points(Points const& points) {
        distance_t x_mean = 0.0, y_mean = 0.0;
        for (Points::size_type i = 1; i <= points.size(); ++i) {
            x_mean += (points[i - 1].x - x_mean) / i;
            y_mean += (points[i - 1].y - y_mean) / i;
        }
        return {x_mean, y_mean};
    }
}

Curve clustering::fsa_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints) {
    std::vector<Points> matchings;
    auto const& center_curve = cluster.center_curve;
    Curve new_center_curve;

    for (auto curve_id: cluster.curve_ids) {
        auto const& curve = curves[curve_id];
        auto matching = frechet::calcMatching(cluster.center_curve, curve);
        matchings.push_back(std::move(matching));
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve[0]);

    for (PointID point_id = fix_endpoints;
            point_id < center_curve.size() - fix_endpoints; ++point_id) {
        Points matching_points;
        for (auto const& matching: matchings) {
            matching_points.push_back(matching[point_id]);
        }
        auto min_enclosing_circle = calcMinEnclosingCircle(matching_points);
        Point new_point = min_enclosing_circle.center;
        if (new_center_curve.empty() ||
                !approx_equal(new_point, new_center_curve.back()))
            new_center_curve.push_back(new_point);
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve.back());
    return new_center_curve;
}

Curve clustering::dba_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints) {
    auto const& center_curve = cluster.center_curve;
    Curve new_center_curve;
    std::vector<std::vector<Points>> matchings;

    for (auto const& curve_id: cluster.curve_ids) {
        auto const& curve = curves[curve_id];
        auto pair_matching = DTW(center_curve, curve).matching();
        std::vector<Points> matching(center_curve.size());

        // Take care of multiple vertices on the curve all mapped to the same
        // vertex of the center.
        std::size_t index = 0;
        for (std::size_t i = 0; i < cluster.center_curve.size(); ++i) {
            while (index < pair_matching.size()
                    && pair_matching[index].first == i) {
                matching[i].emplace_back(curve[pair_matching[index].second]);
                ++index;
            }
        }
        matchings.push_back(std::move(matching));
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve[0]);

    for (PointID i = fix_endpoints;
            i < center_curve.size() - fix_endpoints; ++i) {
        Points points_to_average;

        for (auto const& matching: matchings)
            for (auto const& point: matching[i])
                points_to_average.push_back(point);

        Point new_point = mean_of_points(points_to_average);

        if (new_center_curve.empty() ||
                !approx_equal(new_point, new_center_curve.back()))
            new_center_curve.push_back(new_point);
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve.back());
    return new_center_curve;
}

Curve clustering::cdba_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints) {
    auto const& center_curve = cluster.center_curve;
    Curve new_center_curve;
    std::vector<std::vector<Points>> matchings;

    for (auto const& curve_id: cluster.curve_ids) {
        auto const& curve = curves[curve_id];
        auto const res = std::max(static_cast<std::size_t>(
            (center_curve.curve_length() + curve.curve_length())
            / (center_curve.size() + curve.size()) / 5), 1UL);
        auto param_space_path =  IntegralFrechet(center_curve, curve,
            ParamMetric::L1, res, nullptr).compute_matching().matching;
        std::vector<Points> matching(center_curve.size());

        for (std::size_t i = 0; i < center_curve.size(); ++i) {
            distance_t dist = center_curve.curve_length(i);
            auto y_range = get_y_range(param_space_path, dist);

            if (approx_equal(y_range.first, y_range.second))
                matching[i].emplace_back(curve.interpolate_at(
                    curve.get_cpoint_after(y_range.first)));
            else {
                distance_t length = y_range.second - y_range.first;
                distance_t ratio_to_curve_length = length / curve.curve_length();
                auto number_of_samples = static_cast<std::size_t>(
                    2 * curve.size() * ratio_to_curve_length);
                if (number_of_samples == 0)
                    number_of_samples = 1;

                for (std::size_t j = 0; j < number_of_samples; ++j) {
                    distance_t dist_along_curve = y_range.first
                        + j * length / number_of_samples;
                    matching[i].emplace_back(curve.interpolate_at(
                        curve.get_cpoint_after(dist_along_curve)));
                }
            }
        }
        matchings.push_back(std::move(matching));
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve[0]);

    for (PointID i = fix_endpoints;
            i < center_curve.size() - fix_endpoints; ++i) {
        Points points_to_average;
        for (auto const& matching: matchings)
            for (auto const& point: matching[i])
                points_to_average.emplace_back(point);

        Point new_point = mean_of_points(points_to_average);

        if (new_center_curve.empty() ||
                !approx_equal(new_point, new_center_curve.back()))
            new_center_curve.push_back(new_point);
    }

    if (fix_endpoints)
        new_center_curve.push_back(center_curve.back());
    return new_center_curve;
}

Curve clustering::wedge_update(Curves const& curves, Cluster const& cluster,
        bool fix_endpoints, distance_t eps, int radius) {
    const auto& center_curve = cluster.center_curve;
    Curve new_center_curve;
    std::map<CurveID, Points> matching_paths;

    // We change this point later if fix_endpoints = false.
    new_center_curve.push_back(center_curve[0]);

    for (std::size_t i = 1; i < center_curve.size() - 1; ++i) {
        Points vertices{new_center_curve.back(), center_curve[i],
            center_curve[i + 1]};
        Wedge wedge(vertices, WedgePoints());

        for (auto const& curve_id: cluster.curve_ids) {
            Curve const& curve = curves[curve_id];
            Points param_space_path;
            if (matching_paths.find(curve_id) == matching_paths.end()) {
                auto const res = std::max(static_cast<std::size_t>(
                    (center_curve.curve_length() + curve.curve_length())
                    / (center_curve.size() + curve.size()) / 5), 1UL);
                param_space_path =  IntegralFrechet(center_curve, curve,
                    ParamMetric::L1, res, nullptr).compute_matching().matching;
                matching_paths.emplace(curve_id, param_space_path);
            } else
                param_space_path = matching_paths.at(curve_id);

            WedgePoints seg_1_points = get_points_matched_to_segment(
                param_space_path, center_curve, curve, i - 1, 0);
            WedgePoints seg_2_points = get_points_matched_to_segment(
                param_space_path, center_curve, curve, i, 1);

            wedge.wedge_points.insert(wedge.wedge_points.end(),
                std::make_move_iterator(seg_1_points.begin()),
                std::make_move_iterator(seg_1_points.end()));
            wedge.wedge_points.insert(wedge.wedge_points.end(),
                std::make_move_iterator(seg_2_points.begin()),
                std::make_move_iterator(seg_2_points.end()));
        }

        Point new_point = grid_search(wedge, eps, radius);
        if (!approx_equal(new_point, new_center_curve.back()))
            new_center_curve.push_back(new_point);
    }

    if (fix_endpoints) {
        new_center_curve.push_back(center_curve.back());
        return new_center_curve;
    }

    Points first_points{{0, 0}, new_center_curve[0], new_center_curve[1]};
    Points last_points{new_center_curve.back(), center_curve.back(), {0, 0}};
    Wedge first_wedge(first_points, WedgePoints()),
        last_wedge(last_points, WedgePoints());
    for (auto const& curve_id: cluster.curve_ids) {
        Curve const& curve = curves[curve_id];
        Points param_space_path = matching_paths.at(curve_id);
        WedgePoints first_wps = get_points_matched_to_segment(param_space_path,
            center_curve, curve, 0, 1);
        WedgePoints last_wps = get_points_matched_to_segment(param_space_path,
            center_curve, curve, center_curve.size() - 2, 0);
        first_wedge.wedge_points.insert(first_wedge.wedge_points.end(),
            std::make_move_iterator(first_wps.begin()),
            std::make_move_iterator(first_wps.end()));
        last_wedge.wedge_points.insert(last_wedge.wedge_points.end(),
            std::make_move_iterator(last_wps.begin()),
            std::make_move_iterator(last_wps.end()));
    }
    auto first_point = grid_search(first_wedge, eps, radius);
    auto last_point = grid_search(last_wedge, eps, radius);

    // This is ugly, but what can we do?
    Points final_points(new_center_curve.size());
    final_points[0] = first_point;
    for (std::size_t i = 1; i < new_center_curve.size(); ++i)
        final_points[i] = new_center_curve[i];
    if (!approx_equal(last_point, new_center_curve.back()))
        final_points.push_back(std::move(last_point));
    return Curve(final_points);
}
