#include "Curve.h"
#include <algorithm>
#include <utility>
#include <iostream>

void Curve::set_extreme_points() {
    auto const& front = points.front();
    extreme_points = {front.x, front.y, front.x, front.y};
    prefix_length[0] = 0.0;

    for (PointID i = 1; i < points.size(); ++i) {
        auto const segment_distance = points[i - 1].dist(points[i]);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points.min_x = std::min(extreme_points.min_x, points[i].x);
        extreme_points.min_y = std::min(extreme_points.min_y, points[i].y);
        extreme_points.max_x = std::max(extreme_points.max_x, points[i].x);
        extreme_points.max_y = std::max(extreme_points.max_y, points[i].y);
    }
}

Curve::Curve(std::string name) : m_name(std::move(name)) {}

Curve::Curve(Points pts) : m_name("dummy_name"), points(std::move(pts)),
        prefix_length(points.size()) {
    if (!points.empty())
        set_extreme_points();
}

Curve::Curve(std::string name, Points pts): m_name(std::move(name)),
        points(std::move(pts)), prefix_length(points.size()) {
    if (!points.empty())
        set_extreme_points();
}

void Curve::push_back(Point point) {
    if (prefix_length.empty())
        prefix_length.push_back(0.0);
    else {
        auto const segment_distance = points.back().dist(point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }

    extreme_points.min_x = std::min(extreme_points.min_x, point.x);
    extreme_points.min_y = std::min(extreme_points.min_y, point.y);
    extreme_points.max_x = std::max(extreme_points.max_x, point.x);
    extreme_points.max_y = std::max(extreme_points.max_y, point.y);

    points.push_back(std::move(point));
}

distance_t Curve::curve_length(CPoint const& p) const {
    assert(p.getFraction() >= 0.0 && p.getFraction() <= 1.0);
    assert((p.getPoint() < points.size() - 1)
        || (p.getPoint() == points.size() - 1 && approx_zero(p.getFraction())));
    return approx_zero(p.getFraction())
        ? prefix_length[p.getPoint()]
        : prefix_length[p.getPoint()] * (1. - p.getFraction())
            + prefix_length[p.getPoint() + 1] * p.getFraction();
}

Point Curve::interpolate_at(CPoint const& p) const {
    assert(p.getFraction() >= 0. && p.getFraction() <= 1.);
    assert((p.getPoint() < points.size() - 1)
        || (p.getPoint() == points.size() - 1 && approx_zero(p.getFraction())));
    return approx_zero(p.getFraction())
        ? points[p.getPoint()]
        : points[p.getPoint()] * (1. - p.getFraction())
            + points[p.getPoint() + 1] * p.getFraction();
}

Curve Curve::slice(PointID i, PointID j) const {
    using dtype = decltype(points.begin())::difference_type;
    Points subpoints(points.begin() + static_cast<dtype>(i),
                     points.begin() + static_cast<dtype>(j) + 1);
    return Curve("", subpoints);
}

distance_t Curve::get_fraction(PointID id, distance_t dist) const {
    assert((id < points.size() - 1)
        || (id == points.size() - 1 && approx_zero(dist)));
    assert(approx_zero(dist)
        || (0.0 < dist && dist <= curve_length(id, id + 1) + ABS_TOL));
    return approx_zero(dist) ? 0.0
        : std::min(1.0, dist / curve_length(id, id + 1));
}

CPoint Curve::get_cpoint(PointID id, distance_t dist) const {
    return {id, get_fraction(id, dist)};
}

CPoint Curve::get_cpoint_after(distance_t dist, PointID after_id) const {
    assert(curve_length(after_id) < dist
        || approx_equal(curve_length(after_id), dist));

    while (after_id + 1 < size() && curve_length(after_id + 1) < dist)
        ++after_id;
    return get_cpoint(after_id, dist - curve_length(after_id));
}

SimplifiedCurve Curve::simplify() const {
    std::size_t const k = 10;

    Curve curve(m_name);
    std::vector<PointID> original_points;

    for (PointID i = 0; i < points.size(); ++i) {
        if (i % k == 0 || i + 1 == points.size()) {
            curve.push_back(points[i]);
            original_points.push_back(i);
        }
    }

    return {std::move(curve), std::move(original_points)};
}

Curve Curve::naive_l_simplification(std::size_t l) const {
    Curve curve(m_name);
    std::size_t const period = points.size() > l ? points.size() / l : 1;

    for (std::size_t i = 0; i < points.size(); i += period) {
        if (i == 0 || !approx_equal(curve.back(), points[i]))
            curve.push_back(points[i]);
    }
    if (!approx_equal(curve.back(), points.back()))
        curve.push_back(points.back());
    return curve;
}

auto Curve::getExtremePoints() const -> ExtremePoints const& {
    return extreme_points;
}

distance_t Curve::getUpperBoundDistance(Curve const& other) const {
    auto const& extreme1 = this->getExtremePoints();
    auto const& extreme2 = other.getExtremePoints();

    Point min_point {
        std::min(extreme1.min_x, extreme2.min_x),
        std::min(extreme1.min_y, extreme2.min_y)
    };
    Point max_point {
        std::max(extreme1.max_x, extreme2.max_x),
        std::max(extreme1.max_y, extreme2.max_y)
    };
    return min_point.dist(max_point);
}
