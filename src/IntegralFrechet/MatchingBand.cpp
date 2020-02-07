#include <iostream>
#include "../io.h"
#include "MatchingBand.h"

namespace {
    Point lowest_with_x(const Points& matching, distance_t x, size_t i) {
        // Find index of the first point p with p.x >= x
        while (matching[i].x < x) {
            ++i;

            if (i == matching.size()) {
                return matching[i - 1];
            }
        }

        if (matching[i].x == x) {
            return matching[i];
        }

        distance_t t = (x - matching[i - 1].x) / (matching[i].x - matching[i - 1].x);
        return matching[i - 1] * (1 - t) + matching[i] * t;
    }

    Point lowest_with_y(const Points& matching, distance_t y, size_t i) {
        // Find index of the first point p with p.y >= y
        while (matching[i].y < y) {
            ++i;

            if (i == matching.size()) {
                return matching[i - 1];
            }
        }

        if (matching[i].y == y) {
            return matching[i];
        }

        distance_t t = (y - matching[i - 1].y) / (matching[i].y - matching[i - 1].y);
        return matching[i - 1] * (1 - t) + matching[i] * t;
    }
}

MatchingBand::MatchingBand(const Curve& curve_x, const Curve& curve_y, const Points& matching, distance_t radius) : lower_y_at_x(curve_x.size()), upper_y_at_x(curve_x.size()), lower_x_at_y(curve_y.size()), upper_x_at_y(curve_y.size()) {

    MonotoneComparator compare(BFDirection::Forward);

    // Last observed vertex in matching
    size_t i = 0;

    // Current center of the x/y_interval (= [x/y - radius, x/y + radius])
    Point p(0, 0);

    // First an last points on curve_x within x_interval
    PointID min_x_incl = 0;
    PointID min_x_excl = 0;
    PointID max_x = 0;
    while (max_x + 1 < curve_x.size() && curve_x.curve_length(p.x, max_x + 1) <= radius) {
        max_x += 1;
    }

    // First an last points on curve_y within y_interval
    PointID min_y_incl = 0;
    PointID min_y_excl = 0;
    PointID max_y = 0;
    while (max_y + 1 < curve_y.size() && curve_y.curve_length(p.y, max_y + 1) <= radius) {
        max_y += 1;
    }

    while (!approx_equal(p, matching.back())) {
        // Find next event.
        //
        // Next event is lowest of:
        // - Lowest point in matching with x = min_x_excl (as distance_t) + radius   \ When point is about to leave
        // - Lowest point in matching with y = min_y_excl (as distance_t) + radius   / the interval.
        //
        // - Lowest point in matching with x = (max_x + 1) (as distance_t) - radius     \ When point enters
        // - Lowest point in matching with y = (max_y + 1) (as distance_t) - radius     / the interval.
        //
        // - Next vertex in matching    > When matching takes a turn (separate if-statement, since we need to increment i)
        p = std::min({
            lowest_with_x(matching, curve_x.curve_length(min_x_excl) + radius, i),
            lowest_with_y(matching, curve_y.curve_length(min_y_excl) + radius, i),
        }, compare);

        if (max_x + 1 < curve_x.size()) {
            p = std::min(p, lowest_with_x(matching, curve_x.curve_length(max_x + 1) - radius, i), compare);
        }
        if (max_y + 1 < curve_y.size()) {
            p = std::min(p, lowest_with_y(matching, curve_y.curve_length(max_y + 1) - radius, i), compare);
        }

        if (!compare(p, matching[i + 1])) {
            p = matching[i + 1];
            ++i;
        }

        // Update state
        while (curve_x.curve_length(min_x_incl) < p.x - radius) {
            ++min_x_incl;
        }
        while (curve_y.curve_length(min_y_incl) < p.y - radius) {
            ++min_y_incl;
        }
        while (curve_x.curve_length(min_x_excl) <= p.x - radius + ABS_TOL) {
            ++min_x_excl;
        }
        while (curve_y.curve_length(min_y_excl) <= p.y - radius + ABS_TOL) {
            ++min_y_excl;
        }
        while (max_x + 1 < curve_x.size() && curve_x.curve_length(max_x + 1) <= p.x + radius + ABS_TOL) {
            ++max_x;
        }
        while (max_y + 1 < curve_y.size() && curve_y.curve_length(max_y + 1) <= p.y + radius + ABS_TOL) {
            ++max_y;
        }


        for (PointID ix = min_x_incl; ix <= max_x; ++ix) {
            auto lo = curve_y.get_cpoint(std::clamp(p.y - radius, 0.0, curve_y.curve_length()));
            auto hi = curve_y.get_cpoint(std::clamp(p.y + radius, 0.0, curve_y.curve_length()));

            if (!lower_y_at_x.at(ix).getPoint().valid() || lo < lower_y_at_x.at(ix)) {
                lower_y_at_x.at(ix) = lo;
            }

            if (!upper_y_at_x.at(ix).getPoint().valid() || hi > upper_y_at_x.at(ix)) {
                upper_y_at_x.at(ix) = hi;
            }
        }
        for (PointID iy = min_y_incl; iy <= max_y && iy < curve_y.size(); ++iy) {
            auto lo = curve_x.get_cpoint(std::clamp(p.x - radius, 0.0, curve_x.curve_length()));
            auto hi = curve_x.get_cpoint(std::clamp(p.x + radius, 0.0, curve_x.curve_length()));

            if (!lower_x_at_y.at(iy).getPoint().valid() || lo < lower_x_at_y.at(iy)) {
                lower_x_at_y.at(iy) = lo;
            }

            if (!upper_x_at_y.at(iy).getPoint().valid() || hi > upper_x_at_y.at(iy)) {
                upper_x_at_y.at(iy) = hi;
            }
        }
    }

    Points debug_points;

    for (PointID x = 0; x < curve_x.size(); ++x) {
        debug_points.push_back({
            curve_x.curve_length(x),
            curve_y.curve_length(lower_y_at_x[x]),
        });
        debug_points.push_back({
            curve_x.curve_length(x),
            curve_y.curve_length(upper_y_at_x[x]),
        });
    }

    for (PointID y = 0; y < curve_y.size(); ++y) {
        debug_points.push_back({
            curve_x.curve_length(lower_x_at_y[y]),
            curve_y.curve_length(y),
        });
        debug_points.push_back({
            curve_x.curve_length(upper_x_at_y[y]),
            curve_y.curve_length(y),
        });
    }

    io::export_points("data/out/debug_points.csv", debug_points);
}
