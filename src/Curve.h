#pragma once

#include "Edge.h"

struct SimplifiedCurve;

class Curve {
    const std::string m_name;
    Points points;

    // Total arc length of the curve up to the i-th point
    std::vector<distance_t> prefix_length;

public:
    explicit Curve(std::string name);
    explicit Curve(std::string name, const Points& points);

    std::string name() const {
        return m_name;
    }

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
        return prefix_length.at(id);
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
     * \brief Copy the subcurve c[i], ..., c[j], including both endpoints.
     * \param i The ID of the start of the subcurve.
     * \param j The ID of the last point of the subcurve.
     * \return The copy of the subcurve.
     */
    Curve slice(PointID i, PointID j) const;

    /**
     * Get the fraction that is arc-length dist is along the edge {id, id + 1}.
     */
    distance_t get_fraction(PointID id, distance_t dist) const;

    /**
     * Get the CPoint that is arc-length dist is along the edge {id, id + 1}.
     */
    CPoint get_cpoint(PointID id, distance_t dist) const;

    /**
     * Get the CPoint that is arc-length dist along the curve
     */
    CPoint get_cpoint_after(distance_t dist, PointID after_id = 0) const;

    void push_back(const Point& point);

    Point interpolate_at(const CPoint& point) const;

    SimplifiedCurve simplify() const;
};

struct SimplifiedCurve {
    /**
     * The simplified curve.
     */
    Curve curve;
    /**
     * Map from points in the simplified curve to points in the original curve.
     */
    std::vector<PointID> original_points;
};