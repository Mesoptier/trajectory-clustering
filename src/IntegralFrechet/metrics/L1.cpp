#include "include.h"

namespace {
    void steepest_descent(const Cell& cell, Point s, const Point& t, Points& path) {
        path.push_back(s);

        // If source is monotone greater than target, we need to do steepest
        // descent in reverse monotone direction.
        auto dir = getMonotoneDirection(s, t);
        MonotoneComparator compare(dir);

        // Which return value of Line::side indicates that a point is on the left
        int line_left = (dir == BFDirection::Forward ? -1 : 1);

        auto ell_m_side = cell.ell_m.side(s);

        if (ell_m_side != 0) { // -> s is not on ell_m
            if (ell_m_side == line_left) { // -> s is left of ell_m
                s = std::min({
                    intersect(Line::horizontal(s), cell.ell_m),
                    intersect(Line::horizontal(s), Line::vertical(t))
                }, compare);
            } else  { // -> s is right of ell_m
                s = std::min({
                    intersect(Line::vertical(s), cell.ell_m),
                    intersect(Line::vertical(s), Line::horizontal(t))
                }, compare);
            }

            if (approx_equal(path.back(), s))
                return;
            path.push_back(s);

            ell_m_side = cell.ell_m.side(s);
        }

        if (ell_m_side == 0) { // -> s is on ell_m
            if (compare(s, cell.mid)) {
                // Move along ell_m until center or until the point where steepest descent from t hits ell_m
                s = std::min({
                    cell.mid,
                    intersect(cell.ell_m, Line::vertical(t)),
                    intersect(cell.ell_m, Line::horizontal(t)),
                }, compare);

                if (approx_equal(path.back(), s))
                    return;
                path.push_back(s);
            }
        }
    }
}

template<>
Points compute_matching<ParamMetric::L1>(const Cell& cell) {
    if (cell.is_degenerate()) {
        if (approx_equal(cell.s, cell.t)) {
            return {cell.s};
        } else {
            return {cell.s, cell.t};
        }
    }

    Points path1;
    Points path2;

    path1.reserve(4);
    path2.reserve(4);

    steepest_descent(cell, cell.s, cell.t, path1);
    steepest_descent(cell, cell.t, cell.s, path2);

    // Steepest descent from both ends should find the same minimum
    #ifndef NDEBUG
    if (!approx_equal(path1.back(), path2.back())) {
        throw std::logic_error("paths should find same minimum");
    }
    #endif

    // Combine the two steepest descent paths into one, skipping the duplicated minimum
    path1.insert(path1.end(), path2.rbegin() + 1, path2.rend());
    return path1;
}

template<>
distance_t integrate_linear_dist<ParamMetric::L1>(const Cell& cell) {
    return cell.len1 + cell.len2;
}