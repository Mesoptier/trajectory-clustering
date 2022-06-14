#include "Cell.h"
#include "metrics/include.h"
#include <iostream>

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


distance_t integrand(distance_t mu, std::pair<Point, Point> s, std::pair<Point, Point> t) {

    Point s1 = s.first;
    Point s2 = s.second;
    Point t1 = t.first;
    Point t2 = t.second;

    return abs(s1.x + mu*(t1.x - s1.x) - (s2.x + mu*(t2.x - s2.x))) +
    abs(s1.y + mu*(t1.y - s1.y) - (s2.y + mu*(t2.y - s2.y)));

}

distance_t integrate_L1(const Cell& cell, const Point& s, const Point& t) {
    const auto sp = cell.interpolate_at(s);
    const auto tp = cell.interpolate_at(t);

    distance_t step_count = 1000;
    distance_t interval_length = 1 / step_count;

    distance_t total = 0;

    for (int i = 0; i < step_count; ++i) {
        distance_t mu1 = interval_length * i;
        distance_t mu2 = interval_length * (i + 1);

        total += .5 * interval_length * (integrand(mu2, sp, tp) + integrand(mu1, sp, tp));
    }

    return total;
}


distance_t integrate_linear_cost(const Cell& cell, const Point& s, const Point& t) {
    return integrate_L1(cell, s, t);
    
    const auto[s1, s2] = cell.interpolate_at(s);
    const auto[t1, t2] = cell.interpolate_at(t);

    if (approx_equal(s, Point(0, 0)) && approx_equal(t, Point(1, 1))) {
        int i = 0;
    }

    // Get difference in x- and y-coordinates at the start and end of linear matching (s, t)
    const auto[dx1, dy1] = s1 - s2;
    const auto[dx2, dy2] = t1 - t2;

    const distance_t a = pow(dx1 - dx2, 2) + pow(dy1 - dy2, 2);
    const distance_t b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
    const distance_t c = dx1 * dx1 + dy1 * dy1;

    distance_t val = (a / 3 + b / 2 + c);

    return (a / 3 + b / 2 + c);
}
