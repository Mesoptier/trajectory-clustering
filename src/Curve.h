#ifndef CODE_CURVE_H
#define CODE_CURVE_H

#include <armadillo>
#include "Vertex.h"
#include "Edge.h"

template<class V>
class Curve {

    // Number of dimensions
    unsigned int D;

    // Number of vertices
    unsigned int N;

    // Matrix with dimensions D x N
    arma::Mat<V> vertices;

    // Total arc length of the curve
    V length;

    // Total arc length of the curve up to the i-th vertex
    std::vector<V> lengths;

public:
    explicit Curve(const arma::Mat<V>& vertices)
        : vertices(vertices), D(vertices.n_cols), N(vertices.n_rows), lengths(N) {
        // Compute lengths
        length = 0;

        if (N > 0) {
            arma::Row<V> prevRow;
            arma::Row<V> row = vertices.row(0);

            for (int i = 1; i < N; ++i) {
                prevRow = row;
                row = vertices.row(i);

                // TODO: configure norm
                length += arma::norm(row - prevRow, 2);
                lengths[i] = length;
            }
        }
    }

    arma::Row<V> interpLength(double length) const {
        return interp(length / getLength());
    }

    // Get a point on the curve at normalized distance t
    arma::Row<V> interp(double t) const {
        if (t < 0 || t > 1) {
            throw std::out_of_range("t must be in range [0.0, 1.0]");
        }

        const V targetLength = length * t;

        // Find the first vertex with length greater or equal to the targetLength
        const auto lb = std::lower_bound(lengths.begin(), lengths.end(), targetLength);
        const auto highIndex = lb - lengths.begin();

        if (highIndex == 0) {
            // There is no previous vertex to interpolate with
            return vertices.row(0);
        }

        const V highLength = *lb;
        const auto lowIndex = highIndex - 1;
        const V lowLength = lengths[lowIndex];

        const V highRatio = (targetLength - lowLength) / (highLength - lowLength);
        return vertices.row(lowIndex) * (1 - highRatio) + vertices.row(highIndex) * highRatio;
    }

    unsigned int getNoVertices() const {
        return N;
    }

    V getLength() const {
        return length;
    }

    const arma::Mat<V>& getVertices() const {
        return vertices;
    }

    Vertex<V> getVertex(unsigned int i) const {
        return vertices.row(i);
    }

    Edge<V> getEdge(unsigned int i) const {
        return {getVertex(i), getVertex(i + 1)};
    }
};

#endif //CODE_CURVE_H
