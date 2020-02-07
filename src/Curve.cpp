#include "Curve.h"
#include <algorithm>
#include <utility>

Curve::Curve(std::string name) : m_name(std::move(name)) {}

Curve::Curve(std::string name, const Points& points): m_name(std::move(name)), points(points) {
    if (points.empty()) {
        return;
    }

    prefix_length.reserve(points.size());
    prefix_length.push_back(0);

    for (int i = 1; i < points.size(); ++i) {
        auto segment_distance = points[i - 1].dist(points[i]);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }
}

void Curve::push_back(const Point& point) {
    if (prefix_length.empty()) {
        prefix_length.push_back(0);
    } else {
        auto segment_distance = points.back().dist(point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }

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

distance_t Curve::get_fraction(PointID id, distance_t dist) const {
    assert((id < points.size() - 1) || (id == points.size() - 1 && dist == 0.));
    assert(dist == 0. || (0. < dist && dist <= curve_length(id, id + 1) + ABS_TOL));
    return dist == 0. ? dist : std::min(1., dist / curve_length(id, id + 1));
}

CPoint Curve::get_cpoint(PointID id, distance_t dist) const {
    return {id, get_fraction(id, dist)};
}

CPoint Curve::get_cpoint(distance_t dist) const {
    PointID id = 0;
    while (curve_length(id + 1) < dist) {
        ++id;
    }
    return get_cpoint(id, dist - curve_length(id));
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
