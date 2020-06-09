#include "IntegralFrechet/Cell.h"
#include "IntegralFrechet/metrics/include.h"

std::pair<Point, Point> Cell::interpolate_at(const Point& p) const {
    return {
        ImplicitEdge::interpolate_at(s1, t1, p.x),
        ImplicitEdge::interpolate_at(s2, t2, p.y)
    };
}

Cell Cell::subcell(const Point& start, const Point& end) const {
    const auto[start1, start2] = interpolate_at(start);
    const auto[end1, end2] = interpolate_at(end);
    return {start1, start2, end1, end2};
}

// TODO: Move methods out of Cell namespace
distance_t integrate_linear_cost(const Cell& cell,
        const Point& s, const Point& t) {
    const auto[s1, s2] = cell.interpolate_at(s);
    const auto[t1, t2] = cell.interpolate_at(t);

    // Get difference in x- and y-coordinates at the start and end of linear
    // matching (s, t)
    const auto[dx1, dy1] = s1 - s2;
    const auto[dx2, dy2] = t1 - t2;

    const distance_t a = pow(dx1 - dx2, 2) + pow(dy1 - dy2, 2);
    const distance_t b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    const distance_t c = dx1 * dx1 + dy1 * dy1;

    return (a / 3 + b / 2 + c);
}
