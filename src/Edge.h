#pragma once

#include <utility>
#include <ostream>

#include "util.h"
#include "geom.h"
#include "expressionml.h"

class Edge
{
public:
    const Point first;
    const Point second;

    Point diff;
    distance_t length;

    Line line;

    Edge(const Point& first, const Point& second):
        first(first), second(second), line(Line::fromTwoPoints(first, second)),
        diff(second - first), length(norm(diff)) {}

    /**
     * Get the Point that lies distance `d` along the edge.
     *
     * @param d
     * @return
     */
    Point interpolate_at(distance_t d) const {
//        assert(0 <= d && d <= length);
        return interp(d / length);
    }

    Point interp(distance_t t) const {
        return first * (1 - t) + second * t;
    }

    distance_t param(const Point& point) const {
        // ASSUMPTION: point lies on this edge's infinite line

        if (line.isVertical()) {
            return (first.y - point.y) / (first.y - second.y) * length;
        } else {
            return (first.x - point.x) / (first.x - second.x) * length;
        }
    }

    void writeExpressionML(ExpressionML::Writer& writer) const {
        writer.openFunction("List");
        writer.writePoint(first);
        writer.writePoint(second);
        writer.closeFunction();
    }

    friend std::ostream& operator<<(std::ostream& out, const Edge& edge) {
        out << "V1: " << edge.first << std::endl
            << "V2: " << edge.second << std::endl;
        return out;
    }
};
