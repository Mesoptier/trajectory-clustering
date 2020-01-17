#include "Cell.h"

namespace NewCell {

    std::pair<Point, Point> Cell::interpolate_at(const Point& p) const {
        return {
            ImplicitEdge::interpolate_at(s1, t1, p.x),
            ImplicitEdge::interpolate_at(s2, t2, p.y)
        };
    }

    Cell Cell::subcell(const Point& s, const Point& t) const {
        const auto[s1, s2] = interpolate_at(s);
        const auto[t1, t2] = interpolate_at(t);
        return {s1, s2, t1, t2};
    }

    // TODO: Move methods out of Cell namespace

    distance_t integrate_linear_cost(const Cell& cell) {
        // Get difference in x- and y-coordinates at the start and end of linear matching (s, t)
        const auto[dx1, dy1] = cell.s1 - cell.s2;
        const auto[dx2, dy2] = cell.t1 - cell.t2;

        const distance_t a = pow(dx1 - dx2, 2) + pow(dy1 - dy2, 2);
        const distance_t b = 2 * (dx1 * dx2 - dx1 * dx1 + dy1 * dy2 - dy1 * dy1);
        const distance_t c = dx1 * dx1 + dy1 * dy1;

        return (a / 3 + b / 2 + c);
    }

    //
    // L1
    //

    template<>
    Points compute_matching<ParamMetric::L1>(const Cell& cell) {
        return {};
    }

    template<>
    distance_t integrate_linear_dist<ParamMetric::L1>(const Cell& cell) {
        return cell.len1 + cell.len2;
    }


    //
    // LInfinity (no shortcuts)
    //

    void join_paths(Points& path_forward, const Points& path_backward) {
        if (!approx_equal(path_forward.back(), path_backward.back())) {
            path_forward.push_back(path_backward.back());
        }
        path_forward.insert(path_forward.end(), path_backward.rbegin() + 1, path_backward.rend());
    }

    Point linf_sd_diag(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        const auto line_left = (dir == BFDirection::Forward ? -1 : 1);
        if (cell.ell_v.side(s) != -line_left || cell.ell_h.side(s) != line_left) {
            return s;
        }

        MonotoneComparator compare(dir);
        const Line diag(s, {1, 1});
        auto new_s = std::min({
            intersect(diag, Line::vertical(t)),
            intersect(diag, Line::horizontal(t)),
            intersect(diag, cell.ell_h),
            intersect(diag, cell.ell_v),
        }, compare);

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }

    Point linf_sd_hv(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        const auto line_left = (dir == BFDirection::Forward ? -1 : 1);
        MonotoneComparator compare(dir);

        Point new_s;

        if (cell.ell_v.side(s) != -line_left && cell.ell_h.side(s) == -line_left) {
            if (compare(s, cell.mid)) {
                new_s = std::min({
                    intersect(Line::horizontal(s), cell.ell_v),
                    intersect(Line::horizontal(s), Line::vertical(t)),
                }, compare);
            } else {
                new_s = std::min({
                    intersect(Line::horizontal(s), cell.ell_h),
                    intersect(Line::horizontal(s), Line::vertical(t)),
                }, compare);
            }
        } else if (cell.ell_h.side(s) != line_left && cell.ell_v.side(s) == line_left) {
            if (compare(s, cell.mid)) {
                new_s = std::min({
                    intersect(Line::vertical(s), cell.ell_h),
                    intersect(Line::vertical(s), Line::horizontal(t)),
                }, compare);
            } else {
                new_s = std::min({
                    intersect(Line::vertical(s), cell.ell_v),
                    intersect(Line::vertical(s), Line::horizontal(t)),
                }, compare);
            }
        } else {
            return s;
        }

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }

    Point linf_sd_ell(Points& path, const Cell& cell, const Point& s, const Point& t, BFDirection dir) {
        MonotoneComparator compare(dir);

        if (!compare(s, cell.mid)) {
            return s;
        }

        Point new_s;

        if (cell.ell_v.side(s) == 0) {
            new_s = std::min({
                cell.mid,
                intersect(cell.ell_v, Line::vertical(t)),
                intersect(cell.ell_v, Line::horizontal(t)),
            }, compare);
        } else if (cell.ell_h.side(s) == 0) {
            new_s = std::min({
                cell.mid,
                intersect(cell.ell_h, Line::vertical(t)),
                intersect(cell.ell_h, Line::horizontal(t)),
            }, compare);
        } else {
            return s;
        }

        if (!approx_equal(new_s, s)) {
            path.push_back(new_s);
            return new_s;
        }
        return s;
    }

    template<>
    Points compute_matching<ParamMetric::LInfinity_NoShortcuts>(const Cell& cell) {
        const auto& ell_h = cell.ell_h;
        const auto& ell_v = cell.ell_v;

        // Check if we need to deal with the case where the steepest descent
        // paths would not meet in a single point:
        // TODO: is_opposite_direction essentially checks what the shape of the ellipses is, use that to implement it!
//        if (cell.is_opposite_direction() && !approx_equal(cell.len1, cell.len2)) {
//            if (cell.len2 > cell.len1) {
//                if (ell_h.side(cell.s) != 0 && ell_h.side(cell.t) != 0 && ell_h.side(cell.s) != ell_h.side(cell.t)) {
//                    // Vertical plateau case
//
//                    const Line diag(
//                        cell.s + Point(cell.s.x, cell.s.y + (cell.len2 - cell.len1) / 2),
//                        {1, 1}
//                    );
//
//                    const auto p = intersect(diag, ell_h);
//
//                    const auto path1 = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell.subcell(cell.s, p));
//                    const auto path2 = compute_matching<ParamMetric::LInfinity_NoShortcuts>(cell.subcell(p, cell.t));
//
//
//                    // TODO:
//                    //  1. Compute position of plateau
//                    //  2. Compute steepest descent from s/t to the ends of the plateau
//                    //  3. Combine the paths
//                    return {cell.s, cell.t};
//                }
//            } else {
//                if (ell_v.side(cell.s) != 0 && ell_v.side(cell.t) != 0 && ell_v.side(cell.s) != ell_v.side(cell.t)) {
//                    // Horizontal plateau case
//                    // TODO: Similar to the vertical plateau case above
//                    return {cell.s, cell.t};
//                }
//            }
//        }

        // Otherwise we can simply do steepest descent from both ends
        // TODO:
        //  1. Compute steepest descent from s/t
        //  2. Combine the paths


        // 1. For each endpoint: If we are in a "move diagonally" region
        //      -> Move diagonally until ell_h/ell_v/cell-border
        // 2. If new endpoints have same x and/or y
        //      -> Trivial path, RETURN
        // 3. If not is_opposite_direction
        //      -> Follow ell_h/ell_v, RETURN
        // 4. Resolve corner case, RETURN

        auto s = cell.s;
        auto t = cell.t;

        Points path_forward = {s};
        Points path_backward = {t};

        // Try steepest descent parallel to middle ellipse axis
        s = linf_sd_diag(path_forward, cell, s, t, BFDirection::Forward);
        t = linf_sd_diag(path_backward, cell, t, s, BFDirection::Backward);

        // Return if trivial
        if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
            join_paths(path_forward, path_backward);
            return path_forward;
        }

        // Try steepest descent parallel to cell axes
        s = linf_sd_hv(path_forward, cell, s, t, BFDirection::Forward);
        t = linf_sd_hv(path_backward, cell, t, s, BFDirection::Backward);

        // Return if trivial
        if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
            join_paths(path_forward, path_backward);
            return path_forward;
        }

        // Try steepest descent along other ellipse axes
        s = linf_sd_ell(path_forward, cell, s, t, BFDirection::Forward);
        t = linf_sd_ell(path_backward, cell, t, s, BFDirection::Backward);

        // Return if trivial
        if (approx_equal(s.x, t.x) || approx_equal(s.y, t.y)) {
            join_paths(path_forward, path_backward);
            return path_forward;
        }

        throw std::logic_error("should have found a path by now");
    }

    template<>
    distance_t integrate_linear_dist<ParamMetric::LInfinity_NoShortcuts>(const Cell& cell) {
        return std::max(cell.len1, cell.len2);
    }
}