#include "IntegralFrechet/Cell.h"
#include "IntegralFrechet/metrics/include.h"

std::pair<Point, Point> Cell::interpolate_at(Point const& p) const {
    return {
        ImplicitEdge::interpolate_at(s1, t1, p.x),
        ImplicitEdge::interpolate_at(s2, t2, p.y)
    };
}

Cell Cell::subcell(Point const& start, Point const& end) const {
    auto const [start1, start2] = interpolate_at(start);
    auto const [end1, end2] = interpolate_at(end);
    return {start1, start2, end1, end2};
}

// TODO: Move methods out of Cell namespace
distance_t integrate_linear_cost(Cell const& cell,
        Point const& s, Point const& t) {
    auto const [s1, s2] = cell.interpolate_at(s);
    auto const [t1, t2] = cell.interpolate_at(t);

    // Get difference in x- and y-coordinates at the start and end of linear
    // matching (s, t)
    auto const [dx1, dy1] = s1 - s2;
    auto const [dx2, dy2] = t1 - t2;

    distance_t const a = (dx1 - dx2) * (dx1 - dx2) + (dy1 - dy2) * (dy1 - dy2);
    distance_t const b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    distance_t const c = dx1 * dx1 + dy1 * dy1;

    return (a / 3 + b / 2 + c);
}
