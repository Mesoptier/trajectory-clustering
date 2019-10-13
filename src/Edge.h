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
    V slope;
    V intercept;

    Edge(const Vertex<V>& first, const Vertex<V>& second) : first(first), second(second) {
        diff = second - first;

        if (approx_equal(first(0), second(0))) {
            // Vertical line segment
            slope = INFINITY;
            intercept = first(0);
        } else {
            slope = (second(1) - first(1)) / (second(0) - first(0));
            intercept = first(1) - slope * first(0);
        }
    }

    V getLength() const {
        return arma::norm(first - second, 2);
    }

    Vertex<V> interpLength(double t) const {
        return interp(t / getLength());
    }

    Vertex<V> interp(double t) const {
        return (1 - t) * first + t * second;
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& edge) {
        os << "V1: ";
        edge.first.print(os);
        os << "V2: ";
        edge.second.print(os);
        return os;
    }
};

#endif //CODE_EDGE_H
