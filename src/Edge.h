#ifndef CODE_EDGE_H
#define CODE_EDGE_H

#include <utility>
#include <ostream>
#include <armadillo>
#include "Vertex.h"
#include "util.h"

template<class V>
class Edge {
public:
    const Vertex<V> first;
    const Vertex<V> second;

    arma::Row<V> diff;
    V length;
    V slope;
    V intercept;

    Edge(const Vertex<V>& first, const Vertex<V>& second) : first(first), second(second) {
        diff = second - first;
        length = arma::norm(diff, 2);

        if (approx_equal(first(0), second(0))) {
            // Vertical line segment
            slope = first(1) < second(1) ? INFINITY : -INFINITY;
            intercept = first(0);
        } else {
            slope = (second(1) - first(1)) / (second(0) - first(0));
            intercept = first(1) - slope * first(0);
        }
    }

    Vertex<V> interpLength(V t) const {
        return interp(t / length);
    }

    Vertex<V> interp(V t) const {
        return (1 - t) * first + t * second;
    }

    V param(const Vertex<V>& point) const {
        // ASSUMPTION: point lies on this edge's infinite line

        if (std::isinf(slope)) {
            return (first(1) - point(1)) / (first(1) - second(1));
        } else {
            return (first(0) - point(0)) / (first(0) - second(0));
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& edge) {
        os << "V1: ";
        edge.first.print(os);
        os << "V2: ";
        edge.second.print(os);
        return os;
    }

    friend Vertex<V> intersectInfinite(const Edge& edge1, const Edge& edge2) {
        // ASSUMPTION: edge1 and edge2 are not parallel

        if (std::isinf(edge2.slope)) {
            return intersectInfinite(edge2, edge1);
        }

        V x;
        if (std::isinf(edge1.slope)) {
            x = edge1.intercept;
        } else {
            x = (edge2.intercept - edge1.intercept) / (edge1.slope - edge2.slope);
        }

        return {x, edge2.slope * x + edge2.intercept};
    }
};

#endif //CODE_EDGE_H
