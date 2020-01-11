#pragma once

#include "Edge.h"

class Curve {
    Points points;

    // Total arc length of the curve up to the i-th point
    std::vector<distance_t> prefix_length;

public:
    Curve() = default;
    explicit Curve(const Points& points);

    std::size_t size() const {
        return points.size();
    }
    bool empty() const {
        return points.empty();
    }
    const Point& operator[](PointID id) const { return points[id]; }

    distance_t curve_length() const {
        return prefix_length.back();
    }
    distance_t curve_length(PointID id) const {
        return prefix_length[id];
    }
    distance_t curve_length(PointID i, PointID j) const {
        return curve_length(j) - curve_length(i);
    }
    distance_t curve_length(const CPoint& point) const;
    distance_t curve_length(const CPoint& i, const CPoint j) const {
        return curve_length(j) - curve_length(i);
    }

    Point front() const { return points.front(); }
    Point back() const { return points.back(); }

    const Points& get_points() const {
        return points;
    }

    /**
     * Get the fraction that is arc-length dist is along the edge {id, id + 1}.
     */
    distance_t get_fraction(PointID id, distance_t dist) const;

    /**
     * Get the CPoint that is arc-length dist is along the edge {id, id + 1}.
     */
    CPoint get_cpoint(PointID id, distance_t dist) const;

    void push_back(const Point& point);

    Point interpolate_at(const CPoint& point) const;

    /**
     * Get a coarsened version of the curve.
     * FIXME: Super crappy!
     */
    Curve coarse() const;
};
