#ifndef EDGE_H
#define EDGE_H

#include <utility>
#include <ostream>

#include "utils/util.h"
#include "geom.h"
#include "utils/expressionml.h"

class Edge {
public:
    Point const first;
    Point const second;

    Point diff;
    distance_t length;

    Line line;

    Edge(Point const& f, Point const& s):
        first(f), second(s), diff(s - f), length(norm(diff)),
        line(Line::fromTwoPoints(f, s)) {}

    /**
     * Get the Point that lies distance `d` along the edge.
     */
    Point interpolate_at(distance_t d) const {
        assert(0 <= d && d <= length);
        return interp(d / length);
    }

    Point interp(distance_t t) const {
        return first * (1 - t) + second * t;
    }

    distance_t param(Point const& point) const {
        // ASSUMPTION: point lies on this edge's infinite line
        if (line.isVertical())
            return (first.y - point.y) / (first.y - second.y) * length;
        return (first.x - point.x) / (first.x - second.x) * length;
    }

    void writeExpressionML(ExpressionML::Writer& writer) const {
        writer.openFunction("List");
        writer.writePoint(first);
        writer.writePoint(second);
        writer.closeFunction();
    }

    friend std::ostream& operator<<(std::ostream& out, Edge const& edge) {
        out << "V1: " << edge.first << std::endl
            << "V2: " << edge.second << std::endl;
        return out;
    }
};
#endif
