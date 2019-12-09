#include "Curve.h"

Curve::Curve(const Points& points) : points(points) {
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
    points.push_back(point);

    if (prefix_length.empty()) {
        prefix_length.push_back(0);
    } else {
        auto segment_distance = points.back().dist(point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }
}