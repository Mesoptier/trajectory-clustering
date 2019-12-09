#ifndef CODE_CURVE_H
#define CODE_CURVE_H

#include <armadillo>
#include "Vertex.h"
#include "Edge.h"

class Curve {
    // Matrix with dimensions D x N
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

    distance_t curve_length() const {
        return prefix_length.back();
    }
    distance_t curve_length(int i, int j) const {
        return prefix_length[j] - prefix_length[i];
    }

    Point front() const { return points.front(); }
    Point back() const { return points.back(); }

    void push_back(const Point& point);

    // TODO: Rewrite to match Curve::interpolate_at from klcluster
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
    const Points& getPoints() const {
        return points;
    }

    [[deprecated]]
    Point getPoint(unsigned int i) const {
        return points[i];
    }

    [[deprecated]]
    Edge getEdge(unsigned int i) const {
        return {getPoint(i), getPoint(i + 1)};
    }
};

#endif //CODE_CURVE_H
