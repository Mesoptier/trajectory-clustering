#ifndef CODE_CURVE_H
#define CODE_CURVE_H

#include "Vertex.h"
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

    [[deprecated("Use intepolate_at() instead")]]
    Point interpLength(distance_t length) const {
        // Find the first vertex with length greater or equal to the targetLength
        const auto lb = std::lower_bound(prefix_length.begin(), prefix_length.end(), length);
        const auto highIndex = lb - prefix_length.begin();

        if (highIndex == 0) {
            // There is no previous vertex to interpolate with
            return points[0];
        }

        const distance_t highLength = *lb;
        const auto lowIndex = highIndex - 1;
        const distance_t lowLength = prefix_length[lowIndex];

        const distance_t highRatio = (length - lowLength) / (highLength - lowLength);
        return points[lowIndex] * (1 - highRatio) + points[highIndex] * highRatio;
    }

    [[deprecated("Use curve_length() instead")]]
    distance_t getLength() const {
        return curve_length();
    }

    [[deprecated("Use curve_length(0, i) instead")]]
    distance_t getLength(int i) const {
        return curve_length(0, i);
    }

    [[deprecated]]
    Point getPoint(unsigned int i) const {
        return points[i];
    }

    [[deprecated]]
    Edge getEdge(unsigned int i) const {
        return {getPoint(i), getPoint(i + 1)};
    }

    Edge get_edge(PointID id) const {
        assert(id < points.size() - 1);
        return {points[id], points[id + 1]};
    }
};

#endif //CODE_CURVE_H
