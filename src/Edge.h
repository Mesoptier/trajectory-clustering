#ifndef CODE_EDGE_H
#define CODE_EDGE_H

#include <utility>
#include <ostream>
#include <armadillo>
#include "Vertex.h"
#include "util.h"
#include "geom.h"
#include "expressionml.h"

template<class V>
class Edge {
public:
    const Point first;
    const Point second;

    arma::Row<V> diff;
    V length;

    Line line;

    Edge(const Point& first, const Point& second) : first(first), second(second), line(first, second) {
        diff = second - first;
        length = arma::norm(diff, 2);
    }

    Point interpLength(V t) const {
        return interp(t / length);
    }

    Point interp(V t) const {
        return (1 - t) * first + t * second;
    }

    distance_t param(const Point& point) const {
        // ASSUMPTION: point lies on this edge's infinite line

        if (line.isVertical()) {
            return (first(1) - point(1)) / (first(1) - second(1)) * length;
        } else {
            return (first(0) - point(0)) / (first(0) - second(0)) * length;
        }
    }

    void writeExpressionML(ExpressionML::Writer& writer) const {
        writer.openFunction("List");
        writer.writePoint(first);
        writer.writePoint(second);
        writer.closeFunction();
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& edge) {
        os << "V1: ";
        edge.first.print(os);
        os << "V2: ";
        edge.second.print(os);
        return os;
    }

    [[deprecated("Use intersect on the edge's lines directly instead")]]
    friend Point intersectInfinite(const Edge& edge1, const Edge& edge2) {
        return intersect(edge1.line, edge2.line);
    }
};

#endif //CODE_EDGE_H
