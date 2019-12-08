#ifndef CODE_CURVE_H
#define CODE_CURVE_H

#include <armadillo>
#include "Vertex.h"
#include "Edge.h"

class Curve {
    // Number of vertices
    unsigned int N;

    // Matrix with dimensions D x N
    Points points;

    // Total arc length of the curve
    distance_t length;

    // Total arc length of the curve up to the i-th vertex
    std::vector<distance_t> prefix_length;

public:
    explicit Curve(const Points& points) : points(points), N(points.size()), prefix_length(N) {
        // Compute lengths
        length = 0;

        if (N > 0) {
            Point prevPoint;
            Point point = points[0];

            for (int i = 1; i < N; ++i) {
                prevPoint = point;
                point = points[i];

                length += norm(point - prevPoint);
                prefix_length[i] = length;
            }
        }
    }

    // TODO: Rewrite to match Curve::interpolate_at from klcluster
    Point interpLength(double length) const {
        return interp(length / getLength());
    }

    // Get a point on the curve at normalized distance t
    Point interp(double t) const {
        if (t < 0 || t > 1) {
            throw std::out_of_range("t must be in range [0.0, 1.0]");
        }

        const distance_t targetLength = length * t;

        // Find the first vertex with length greater or equal to the targetLength
        const auto lb = std::lower_bound(prefix_length.begin(), prefix_length.end(), targetLength);
        const auto highIndex = lb - prefix_length.begin();

        if (highIndex == 0) {
            // There is no previous vertex to interpolate with
            return points[0];
        }

        const distance_t highLength = *lb;
        const auto lowIndex = highIndex - 1;
        const distance_t lowLength = prefix_length[lowIndex];

        const distance_t highRatio = (targetLength - lowLength) / (highLength - lowLength);
        return points[lowIndex] * (1 - highRatio) + points[highIndex] * highRatio;
    }

    unsigned int getNoVertices() const {
        return N;
    }

    distance_t getLength() const {
        return length;
    }

    distance_t getLength(int i) const {
        return prefix_length[i];
    }

    const Points& getPoints() const {
        return points;
    }

    Point getPoint(unsigned int i) const {
        return points[i];
    }

    Edge getEdge(unsigned int i) const {
        return {getPoint(i), getPoint(i + 1)};
    }
};

#endif //CODE_CURVE_H
