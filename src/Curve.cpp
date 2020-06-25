#include "Curve.h"
#include <algorithm>
#include <utility>
#include <iostream>

Curve::Curve(std::string name) : m_name(std::move(name)) {}

Curve::Curve(const Points& pts) : m_name("name"), points(pts), prefix_length(pts.size()) {

    if (points.empty()) { return; }

	auto const& front = points.front();
	extreme_points = { front.x, front.y, front.x, front.y };
	prefix_length[0] = 0;


    for (PointID i = 1; i < points.size(); ++i)
	{
		auto segment_distance = points[i - 1].dist(points[i]);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

		extreme_points.min_x = std::min(extreme_points.min_x, points[i].x);
		extreme_points.min_y = std::min(extreme_points.min_y, points[i].y);
		extreme_points.max_x = std::max(extreme_points.max_x, points[i].x);
		extreme_points.max_y = std::max(extreme_points.max_y, points[i].y);
	}
}

Curve::Curve(std::string name, const Points& pts): m_name(std::move(name)), points(pts) {
    if (points.empty()) {
        return;
    }

    prefix_length.reserve(points.size());
    prefix_length.push_back(0);

    for (size_t i = 1; i < points.size(); ++i) {
        auto segment_distance = points[i - 1].dist(points[i]);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }

    for (PointID i = 1; i < points.size(); ++i)
	{
		auto segment_distance = points[i - 1].dist(points[i]);
		prefix_length[i] = prefix_length[i - 1] + segment_distance;

		extreme_points.min_x = std::min(extreme_points.min_x, points[i].x);
		extreme_points.min_y = std::min(extreme_points.min_y, points[i].y);
		extreme_points.max_x = std::max(extreme_points.max_x, points[i].x);
		extreme_points.max_y = std::max(extreme_points.max_y, points[i].y);
	}
}

void Curve::push_back(const Point& point) {
    if (prefix_length.empty()) {
        prefix_length.push_back(0);
    } else {
        auto segment_distance = points.back().dist(point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }

    extreme_points.min_x = std::min(extreme_points.min_x, point.x);
	extreme_points.min_y = std::min(extreme_points.min_y, point.y);
	extreme_points.max_x = std::max(extreme_points.max_x, point.x);
	extreme_points.max_y = std::max(extreme_points.max_y, point.y);

    points.push_back(point);
}

distance_t Curve::curve_length(const CPoint& p) const {
    assert(p.getFraction() >= 0. && p.getFraction() <= 1.);
    assert((p.getPoint() < points.size() - 1) || (p.getPoint() == points.size() - 1 && p.getFraction() == 0.));
    return p.getFraction() == 0.
        ? prefix_length[p.getPoint()]
        : prefix_length[p.getPoint()] * (1. - p.getFraction()) + prefix_length[p.getPoint() + 1] * p.getFraction();
}

Point Curve::interpolate_at(const CPoint& p) const {
    assert(p.getFraction() >= 0. && p.getFraction() <= 1.);
    assert((p.getPoint() < points.size() - 1) || (p.getPoint() == points.size() - 1 && p.getFraction() == 0.));
    return p.getFraction() == 0.
        ? points[p.getPoint()]
        : points[p.getPoint()] * (1. - p.getFraction()) + points[p.getPoint() + 1] * p.getFraction();
}

Curve Curve::slice(PointID i, PointID j) const {
    using dtype = decltype(points.begin())::difference_type;
    Points subpoints(points.begin() + static_cast<dtype>(i),
                     points.begin() + static_cast<dtype>(j) + 1);
    return Curve("", subpoints);
}

distance_t Curve::get_fraction(PointID id, distance_t dist) const {
    // if (dist < 1e-10) {
    //     dist = 0.;
    // }
    assert((id < points.size() - 1) || (id == points.size() - 1 && dist == 0.));
    assert(dist == 0. || (0. < dist && dist <= curve_length(id, id + 1) + ABS_TOL));
    return dist == 0. ? dist : std::min(1., dist / curve_length(id, id + 1));
}

CPoint Curve::get_cpoint(PointID id, distance_t dist) const {
    return {id, get_fraction(id, dist)};
}

CPoint Curve::get_cpoint_after(distance_t dist, PointID after_id) const {
    assert(curve_length(after_id) <= dist);

    while (after_id + 1 < size() && curve_length(after_id + 1) < dist) {
        ++after_id;
    }
    
    return get_cpoint(after_id, dist - curve_length(after_id));
}

Curve Curve::simplify(bool maintain_lengths) const {
    Curve other(m_name);

    for (size_t i = 0; i < points.size(); ++i) {
        if (i % 10 == 0 || i == points.size() - 1) {
            if (maintain_lengths) {
                other.points.push_back(points[i]);
                other.prefix_length.push_back(prefix_length[i]);
            } else {
                other.push_back(points[i]);
            }
        }
    }

    return other;
}

Curve Curve::naive_l_simplification(int l) const {
    Curve other(m_name);

    size_t period;

    if (points.size() > l) {
        period = points.size() / l;
    } else {
        period = 1;
    }

    for (size_t i = 0; i < points.size(); i += period) {
        if (other.get_points().empty())
            other.push_back(points[i]);
        else if(!approx_equal(other.get_points().back(), points[i])) {
            other.push_back(points[i]);
        }
    }

    if (other.get_points().back() != points.back()) {
        other.push_back(points.back());
    }

    return other;
}

auto Curve::getExtremePoints() const -> ExtremePoints const&
{
	return extreme_points;
}

distance_t Curve::getUpperBoundDistance(Curve const& other) const
{
	auto const& extreme1 = this->getExtremePoints();
	auto const& extreme2 = other.getExtremePoints();

	Point min_point{ std::min(extreme1.min_x, extreme2.min_x),
		std::min(extreme1.min_y, extreme2.min_y) };
	Point max_point = { std::max(extreme1.max_x, extreme2.max_x),
		std::max(extreme1.max_y, extreme2.max_y) };

	return min_point.dist(max_point);
}
